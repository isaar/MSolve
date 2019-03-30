using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.Commons;

namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// Free dofs are assigned global / subdomain indices in a node major fashion: The dofs of the first node are numbered, then 
    /// the dofs of the second node, etc. Note that the dofs of each node are assumed to be the same and supplied by the client. 
    /// Based on that assumption, this class is much faster than its alternatives. Constrained dofs are ignored.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class UniformDofOrderingStrategy : IFreeDofOrderingStrategy
    {
        private readonly IReadOnlyList<DOFType> dofsPerNode;

        public UniformDofOrderingStrategy(IReadOnlyList<DOFType> dofsPerNode)
        {
            this.dofsPerNode = dofsPerNode;
        }

        public (int numGlobalFreeDofs, DofTable globalFreeDofs) OrderGlobalDofs(IStructuralModel_v2 model)
            => OrderFreeDofsOfNodeSet(model.Nodes, model.Constraints);


        public (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) OrderSubdomainDofs(ISubdomain_v2 subdomain)
            => OrderFreeDofsOfNodeSet(subdomain.Nodes, subdomain.Constraints);


        private (int numFreeDofs, DofTable freeDofs) OrderFreeDofsOfNodeSet(IEnumerable<INode> sortedNodes,
            Table<INode, DOFType, double> constraints)
        {
            var freeDofs = new DofTable();
            int dofCounter = 0;
            foreach (INode node in sortedNodes)
            {
                bool isNodeConstrained = constraints.TryGetDataOfRow(node,
                    out IReadOnlyDictionary<DOFType, double> constraintsOfNode);
                foreach (DOFType dof in dofsPerNode)
                {
                    bool isDofConstrained = isNodeConstrained ? constraintsOfNode.ContainsKey(dof) : false;
                    if (!isDofConstrained) freeDofs[node, dof] = dofCounter++;
                }
            }
            return (dofCounter, freeDofs);
        }
    }
}
