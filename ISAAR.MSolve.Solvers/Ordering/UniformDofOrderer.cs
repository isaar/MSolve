using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.Commons;

namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// Free dofs are assigned global (actually subdomain) indices in a node major fashion: The dofs of the first node are 
    /// numbered, then the dofs of the second node, etc. Contrary to <see cref="NodeMajorDofOrderer"/>, the dofs of each node
    /// are assumed to be the same and supplied by the client. Based on that assumption, this class is much faster than
    /// <see cref="NodeMajorDofOrderer"/>. Constrained dofs are ignored.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class UniformDofOrderer: IDofOrderer
    {
        private readonly IReadOnlyList<DOFType> dofsPerNode;

        public UniformDofOrderer(IReadOnlyList<DOFType> dofsPerNode)
        {
            this.dofsPerNode = dofsPerNode;
        }

        public IGlobalFreeDofOrdering OrderDofs(IStructuralModel_v2 model)
        {
            //TODO: move this to the end
            (int numGlobalFreeDofs, DofTable globalFreeDofs) =
                   OrderFreeDofsOfNodeSet(model.Nodes, model.Constraints);

            // Order subdomain dofs
            var subdomainOrderings = new Dictionary<ISubdomain_v2, ISubdomainFreeDofOrdering>(model.Subdomains.Count);
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) =
                    OrderFreeDofsOfNodeSet(subdomain.Nodes, subdomain.Constraints);
                ISubdomainFreeDofOrdering subdomainOrdering =
                    new SubdomainFreeDofOrderingGeneral(numSubdomainFreeDofs, subdomainFreeDofs, globalFreeDofs);
                subdomainOrderings.Add(subdomain, subdomainOrdering);
            }

            // Order global dofs
            return new GlobalFreeDofOrderingGeneral(numGlobalFreeDofs, globalFreeDofs, subdomainOrderings);
        }

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
