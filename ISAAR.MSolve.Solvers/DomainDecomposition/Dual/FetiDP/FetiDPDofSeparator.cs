using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;

//TODO: Remove code duplication between this and Feti1DofSeparator
//TODO: Perhaps I should also find and expose the indices of boundary remainder and internal remainder dofs into the sequence 
//      of all free dofs of each subdomain
//TODO: Decide which of these data structures will be stored and which will be used ONCE to create all required mapping matrices.
//TODO: Perhaps the corner dof logic should be moved to another class.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP
{
    public class FetiDPDofSeparator : DofSeparatorBase
    {
        /// <summary>
        /// Indices of boundary remainder dofs into the sequence of all remainder dofs of each subdomain.
        /// </summary>
        public override Dictionary<int, int[]> BoundaryDofIndices { get; protected set; }

        /// <summary>
        /// (Node, IDofType) pairs for each boundary remainder dof of each subdomain. Their order is the same one as 
        /// <see cref="BoundaryDofIndices"/>.
        /// </summary>
        public override Dictionary<int, (INode node, IDofType dofType)[]> BoundaryDofs { get; protected set; }

        /// <summary>
        /// Also called Bc in papers by Farhat or Lc in NTUA theses. 
        /// </summary>
        public Dictionary<int, UnsignedBooleanMatrix> CornerBooleanMatrices { get; private set; } //TODO: This should be sparse

        /// <summary>
        /// Indices of (boundary) corner dofs into the sequence of all free dofs of each subdomain.
        /// </summary>
        public Dictionary<int, int[]> CornerDofIndices { get; private set; }

        /// <summary>
        /// Dof ordering for corner dofs of the model: Each (INode, IDofType) pair is associated with the index of that dof into 
        /// a vector corresponding to all corner dofs of the model.
        /// </summary>
        public DofTable GlobalCornerDofOrdering { get; private set; }

        /// <summary>
        /// If Xf is a vector with all free dofs of the model and Xc is a vector with all corner dofs of the model, then
        /// Xf[GlobalCornerToFreeDofMap[i]] = Xc[i].
        /// </summary>
        public int[] GlobalCornerToFreeDofMap { get; set; }

        /// <summary>
        /// Indices of internal remainder dofs into the sequence of all remainder dofs of each subdomain.
        /// </summary>
        public override Dictionary<int, int[]> InternalDofIndices { get; protected set; }

        /// <summary>
        /// The number of corner dofs of the model.
        /// </summary>
        public int NumGlobalCornerDofs { get; private set; }

        /// <summary>
        /// Dof ordering for remainder (boundary and internal) dofs of each subdomain: Each (INode, IDofType) pair of a 
        /// subdomain is associated with the index of that dof into a vector corresponding to remainder dofs of that subdomain.
        /// </summary>
        public Dictionary<int, DofTable> RemainderDofOrderings { get; private set; }

        /// <summary>
        /// Indices of remainder (boundary and internal) dofs into the sequence of all free dofs of each subdomain.
        /// </summary>
        public Dictionary<int, int[]> RemainderDofIndices { get; private set; }

        /// <summary>
        /// Dof ordering for corner dofs of each subdomain: Each (INode, IDofType) pair of a subdomain is associated with the   
        /// index of that dof into a vector corresponding to corner dofs of that subdomain.
        /// </summary>
        public Dictionary<int, DofTable> SubdomainCornerDofOrderings { get; private set; }

        public FetiDPDofSeparator()
        {
            CornerDofIndices = new Dictionary<int, int[]>();
            RemainderDofIndices = new Dictionary<int, int[]>();
            InternalDofIndices = new Dictionary<int, int[]>();
            BoundaryDofIndices = new Dictionary<int, int[]>();
            BoundaryDofs = new Dictionary<int, (INode node, IDofType dofType)[]>();
            RemainderDofOrderings = new Dictionary<int, DofTable>();
            SubdomainCornerDofOrderings = new Dictionary<int, DofTable>();
            CornerBooleanMatrices = new Dictionary<int, UnsignedBooleanMatrix>();
        }

        /// <summary>
        /// Bc unsigned boolean matrices that map global to subdomain corner dofs. This method must be called after 
        /// <see cref="DefineGlobalCornerDofs(IStructuralModel, Dictionary{int, HashSet{INode}})"/>.
        /// </summary>
        public void CalcCornerMappingMatrices(IStructuralModel model, Dictionary<int, HashSet<INode>> subdomainCornerNodes)
        { //TODO: Can I reuse subdomain data? Yes if the global corner dofs have not changed.
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int s = subdomain.ID;
                DofTable localCornerDofOrdering = SubdomainCornerDofOrderings[s];
                int numLocalCornerDofs = localCornerDofOrdering.EntryCount;
                var Bc = new UnsignedBooleanMatrix(numLocalCornerDofs, NumGlobalCornerDofs);
                foreach ((INode node, IDofType dofType, int localIdx) in localCornerDofOrdering)
                {
                    int globalIdx = GlobalCornerDofOrdering[node, dofType];
                    Bc.AddEntry(localIdx, globalIdx);
                }
                CornerBooleanMatrices[s] = Bc;
            }
        }

        public void DefineGlobalBoundaryDofs(IStructuralModel model, Dictionary<int, HashSet<INode>> subdomainCornerNodes)
        {
            //TODO: These might be needed elsewhere too, in which case it should probably be sorted.
            var globalCornerNodes = new HashSet<INode>();
            foreach (IEnumerable<INode> subdomainNodes in subdomainCornerNodes.Values)
            {
                globalCornerNodes.UnionWith(subdomainNodes);
            }
            IEnumerable<INode> globalRemainderNodes = model.Nodes.Where(node => !globalCornerNodes.Contains(node));
            base.DefineGlobalBoundaryDofs(globalRemainderNodes, model.GlobalDofOrdering); //TODO: This could also be reused in some cases
        }

        public void DefineGlobalCornerDofs(IStructuralModel model, Dictionary<int, HashSet<INode>> subdomainCornerNodes)
        {
            // Gather all corner nodes
            //TODO: This is also calculated in SeparateDofs(). Reuse it.
            var globalCornerNodes = new SortedSet<INode>(); //TODO: Can this be optimized?
            foreach (IEnumerable<INode> subdomainNodes in subdomainCornerNodes.Values)
            {
                foreach (INode node in subdomainNodes) globalCornerNodes.Add(node);
            }

            // Order global corner dofs and create the global corner to global free map.
            var cornerToGlobalDofs = new List<int>(globalCornerNodes.Count * 3);
            var globalCornerDofOrdering = new DofTable(); //TODO: Should this be cached?
            int cornerDofCounter = 0;
            foreach (INode cornerNode in globalCornerNodes)
            {
                bool hasFreeDofs = model.GlobalDofOrdering.GlobalFreeDofs.TryGetDataOfRow(cornerNode,
                    out IReadOnlyDictionary<IDofType, int> dofsOfNode);
                if (!hasFreeDofs) throw new Exception($"Corner node {cornerNode.ID} has only constrained or embedded dofs.");
                foreach (var dofTypeIdxPair in dofsOfNode)
                {
                    IDofType dofType = dofTypeIdxPair.Key;
                    int globalDofIdx = dofTypeIdxPair.Value;
                    globalCornerDofOrdering[cornerNode, dofType] = cornerDofCounter++;
                    cornerToGlobalDofs.Add(globalDofIdx);
                }
            }
            NumGlobalCornerDofs = cornerDofCounter;
            GlobalCornerToFreeDofMap = cornerToGlobalDofs.ToArray();
            GlobalCornerDofOrdering = globalCornerDofOrdering;
        }

        /// <summary>
        /// This must be called after <see cref="SeparateCornerRemainderDofs(ISubdomain, HashSet{INode}, IEnumerable{INode})"/>.
        /// </summary>
        public void SeparateBoundaryInternalDofs(ISubdomain subdomain, IEnumerable<INode> remainderAndConstrainedNodes)
        {
            int s = subdomain.ID;

            (int[] internalDofIndices, int[] boundaryDofIndices, (INode node, IDofType dofType)[] boundaryDofConnectivities)
                = DofSeparatorBase.SeparateBoundaryInternalDofs(remainderAndConstrainedNodes, RemainderDofOrderings[s]);
            InternalDofIndices[s] = internalDofIndices;
            BoundaryDofIndices[s] = boundaryDofIndices;
            BoundaryDofs[s] = boundaryDofConnectivities;
        }

        public void SeparateCornerRemainderDofs(ISubdomain subdomain, HashSet<INode> cornerNodes, 
            IEnumerable<INode> remainderAndConstrainedNodes)
        {
            int s = subdomain.ID;

            var cornerDofs = new List<int>();
            var remainderDofs = new List<int>();
            foreach (INode node in cornerNodes)
            {
                IEnumerable<int> dofsOfNode = subdomain.FreeDofOrdering.FreeDofs.GetValuesOfRow(node);
                cornerDofs.AddRange(dofsOfNode);
            }
            foreach (INode node in remainderAndConstrainedNodes)
            {
                IEnumerable<int> dofsOfNode = subdomain.FreeDofOrdering.FreeDofs.GetValuesOfRow(node);
                remainderDofs.AddRange(dofsOfNode);
            }
            CornerDofIndices[s] = cornerDofs.ToArray();
            RemainderDofIndices[s] = remainderDofs.ToArray();

            //TODO: This dof ordering should be optimized, such that the factorization of Krr is efficient.
            RemainderDofOrderings[s] = subdomain.FreeDofOrdering.FreeDofs.GetSubtableForNodes(remainderAndConstrainedNodes);
            SubdomainCornerDofOrderings[s] = subdomain.FreeDofOrdering.FreeDofs.GetSubtableForNodes(cornerNodes);
        }
    }
}
