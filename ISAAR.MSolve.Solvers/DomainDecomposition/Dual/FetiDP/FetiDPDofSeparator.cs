using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: Remove code duplication between this and Feti1DofSeparator
//TODO: Perhaps I should also find and expose the indices of boundary remainder and internal remainder dofs into the sequence 
//      of all free dofs of each subdomain
//TODO: Decide which of these data structures will be cached and which will be used once to create all required mapping matrices.
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
        /// Indices of (boundary) corner dofs into the sequence of all free dofs of each subdomain.
        /// </summary>
        public Dictionary<int, int[]> CornerDofIndices { get; private set; }

        /// <summary>
        /// Dof ordering for corner dofs of each subdomain.
        /// </summary>
        public Dictionary<int, DofTable> SubdomainCornerDofOrderings { get; private set; }

        /// <summary>
        /// Also called Bc in papers by Farhat. 
        /// </summary>
        public Dictionary<int, Matrix> BooleanCornerMatrices { get; private set; } //TODO: This should be sparse

        /// <summary>
        /// Indices of internal remainder dofs into the sequence of all remainder dofs of each subdomain.
        /// </summary>
        public override Dictionary<int, int[]> InternalDofIndices { get; protected set; }

        public int NumGlobalCornerDofs { get; private set; }

        /// <summary>
        /// Dof ordering for remainder (boundary and internal) dofs of each subdomain.
        /// </summary>
        public Dictionary<int, DofTable> RemainderDofOrderings { get; private set; }

        /// <summary>
        /// Indices of remainder (boundary and internal) dofs into the sequence of all free dofs of each subdomain.
        /// </summary>
        public Dictionary<int, int[]> RemainderDofIndices { get; private set; }

        public void DefineCornerMappingMatrices(IStructuralModel model, Dictionary<int, INode[]> subdomainCornerNodes)
        {
            // Gather all corner nodes
            //TODO: This is also calculated in SeparateDofs(). Reuse it.
            var globalCornerNodes = new SortedSet<INode>(); //TODO: Can this be optimized?
            foreach (IReadOnlyList<INode> subdomainNodes in subdomainCornerNodes.Values)
            {
                foreach (INode node in subdomainNodes) globalCornerNodes.Add(node);
            }

            // Order global corner dofs
            var globalCornerDofOrdering = new DofTable(); //TODO: Should this be cached?
            int cornerDofCounter = 0;
            foreach (INode cornerNode in globalCornerNodes)
            {
                IEnumerable<IDofType> dofsOfNode = model.GlobalDofOrdering.GlobalFreeDofs.GetColumnsOfRow(cornerNode);
                Debug.Assert(dofsOfNode.Count() >= 1);
                foreach (IDofType dofType in dofsOfNode)
                {
                    globalCornerDofOrdering[cornerNode, dofType] = cornerDofCounter++;
                }
            }
            NumGlobalCornerDofs = cornerDofCounter;

            // Order local corner dofs
            SubdomainCornerDofOrderings = new Dictionary<int, DofTable>();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                DofTable localCornerDofOrdering =
                    subdomain.FreeDofOrdering.FreeDofs.GetSubtableForNodes(subdomainCornerNodes[subdomain.ID]);
                SubdomainCornerDofOrderings[subdomain.ID] = localCornerDofOrdering;
            }

            // Fill Bc matrix of each subdomain 
            BooleanCornerMatrices = new Dictionary<int, Matrix>();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                DofTable localCornerDofOrdering = SubdomainCornerDofOrderings[subdomain.ID];
                int numLocalCornerDofs = localCornerDofOrdering.EntryCount;
                var Bc = Matrix.CreateZero(numLocalCornerDofs, NumGlobalCornerDofs);
                foreach ((INode node, IDofType dofType, int localIdx) in localCornerDofOrdering)
                {
                    int globalIdx = globalCornerDofOrdering[node, dofType];
                    Bc[localIdx, globalIdx] = 1;
                }
                BooleanCornerMatrices[subdomain.ID] = Bc;
            }
        }

        public void SeparateDofs(IStructuralModel model, Dictionary<int, INode[]> subdomainCornerNodes)
        {
            //TODO: These might be needed elsewhere too, in which case it should probably be sorted.
            var allCornerNodes = new HashSet<INode>();
            foreach (IReadOnlyList<INode> subdomainNodes in subdomainCornerNodes.Values)
            {
                foreach (INode node in subdomainNodes) allCornerNodes.Add(node);
            }
            IEnumerable<INode> allRemainderNodes = model.Nodes.Where(node => !allCornerNodes.Contains(node));

            base.GatherGlobalBoundaryDofs(allRemainderNodes, model.GlobalDofOrdering);

            CornerDofIndices = new Dictionary<int, int[]>();
            RemainderDofIndices = new Dictionary<int, int[]>();
            InternalDofIndices = new Dictionary<int, int[]>();
            BoundaryDofIndices = new Dictionary<int, int[]>();
            BoundaryDofs = new Dictionary<int, (INode node, IDofType dofType)[]>();
            RemainderDofOrderings = new Dictionary<int, DofTable>();

            foreach (ISubdomain subdomain in model.Subdomains)
            {
                var cornerNodes = new HashSet<INode>(subdomainCornerNodes[subdomain.ID]);
                INode[] remainderAndConstrainedNodes = subdomain.Nodes.Where(node => !cornerNodes.Contains(node)).ToArray(); //TODO: extract this

                // Separate corner / remainder dofs
                var cornerDofs = new List<int>();
                var remainderDofs = new List<int>();
                foreach (INode node in subdomainCornerNodes[subdomain.ID])
                {
                    IEnumerable<int> dofsOfNode = subdomain.FreeDofOrdering.FreeDofs.GetValuesOfRow(node);
                    cornerDofs.AddRange(dofsOfNode);
                }
                foreach (INode node in remainderAndConstrainedNodes)
                {
                    IEnumerable<int> dofsOfNode = subdomain.FreeDofOrdering.FreeDofs.GetValuesOfRow(node);
                    remainderDofs.AddRange(dofsOfNode);
                }
                CornerDofIndices[subdomain.ID] = cornerDofs.ToArray();
                RemainderDofIndices[subdomain.ID] = remainderDofs.ToArray();

                //TODO: These dof ordering should be optimized, such that the factorization of Krr is efficient.
                DofTable remainderDofOrdering = subdomain.FreeDofOrdering.FreeDofs.GetSubtableForNodes(remainderAndConstrainedNodes);

                // Separate internal / boundary dofs
                (int[] internalDofIndices, int[] boundaryDofIndices, (INode node, IDofType dofType)[] boundaryDofConnectivities)
                    = base.SeparateBoundaryInternalDofs(remainderAndConstrainedNodes, remainderDofOrdering);

                InternalDofIndices[subdomain.ID] = internalDofIndices;
                BoundaryDofIndices[subdomain.ID] = boundaryDofIndices;
                BoundaryDofs[subdomain.ID] = boundaryDofConnectivities;
                RemainderDofOrderings[subdomain.ID] = remainderDofOrdering;
            }
        }
    }
}
