using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: Remove code duplication between this and Feti1DofSeparator
//TODO: Perhaps I should also find and expose the indices of boundary remainder and internal remainder dofs into the sequence 
//      of all free dofs of each subdomain
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP
{
    public class FetiDPDofSeparator : DofSeparatorBase
    {
        internal Dictionary<int, (INode node, IDofType dofType)[]> BoundaryDofConnectivities { get; private set; }
        internal Dictionary<int, int[]> BoundaryDofMultiplicities { get; private set; }

        /// <summary>
        /// Indices of boundary remainder dofs into the sequence of all remainder dofs of each subdomain.
        /// </summary>
        internal Dictionary<int, int[]> BoundaryIntoRemainderDofIndices { get; private set; }

        /// <summary>
        /// Indices of (boundary) corner dofs into the sequence of all free dofs of each subdomain.
        /// </summary>
        internal Dictionary<int, int[]> CornerIntoFreeDofIndices { get; private set; }

        /// <summary>
        /// Indices of internal remainder dofs into the sequence of all remainder dofs of each subdomain.
        /// </summary>
        internal Dictionary<int, int[]> InternalIntoRemainderDofIndices { get; private set; }

        /// <summary>
        /// Indices of remainder (boundary and internal) dofs into the sequence of all free dofs of each subdomain.
        /// </summary>
        internal Dictionary<int, int[]> RemainderIntoFreeDofIndices { get; private set; }

        internal void SeparateDofs(IStructuralModel model, Dictionary<int, INode[]> subdomainCornerNodes)
        {
            //TODO: These might be needed elsewhere too, in which case it should probably be sorted.
            var allCornerNodes = new HashSet<INode>();
            foreach (IReadOnlyList<INode> subdomainNodes in subdomainCornerNodes.Values)
            {
                foreach (INode node in subdomainNodes) allCornerNodes.Add(node);
            }
            IEnumerable<INode> allRemainderNodes = model.Nodes.Where(node => !allCornerNodes.Contains(node));

            base.GatherDualDofs(allRemainderNodes, model.GlobalDofOrdering);

            CornerIntoFreeDofIndices = new Dictionary<int, int[]>();
            RemainderIntoFreeDofIndices = new Dictionary<int, int[]>();
            InternalIntoRemainderDofIndices = new Dictionary<int, int[]>();
            BoundaryIntoRemainderDofIndices = new Dictionary<int, int[]>();
            BoundaryDofMultiplicities = new Dictionary<int, int[]>();
            BoundaryDofConnectivities = new Dictionary<int, (INode node, IDofType dofType)[]>();

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
                CornerIntoFreeDofIndices[subdomain.ID] = cornerDofs.ToArray();
                RemainderIntoFreeDofIndices[subdomain.ID] = remainderDofs.ToArray();

                // Separate internal / boundary dofs

                //TODO: This dof ordering should be optimized, such that the factorization of Krr is efficient. It should also be
                //      cached and reused later.
                DofTable remainderDofsTable = subdomain.FreeDofOrdering.FreeDofs.GetSubtableForNodes(remainderAndConstrainedNodes);

                (int[] internalDofIndices, int[] boundaryDofIndices, int[] boundaryDofMultiplicities,
                    (INode node, IDofType dofType)[] boundaryDofConnectivities) =
                    base.SeparateBoundaryInternalDofs(remainderAndConstrainedNodes, remainderDofsTable);

                InternalIntoRemainderDofIndices[subdomain.ID] = internalDofIndices;
                BoundaryIntoRemainderDofIndices[subdomain.ID] = boundaryDofIndices;
                BoundaryDofMultiplicities[subdomain.ID] = boundaryDofMultiplicities;
                BoundaryDofConnectivities[subdomain.ID] = boundaryDofConnectivities;
            }
        }
    }
}
