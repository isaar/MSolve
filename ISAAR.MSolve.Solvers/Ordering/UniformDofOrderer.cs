using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// Free dofs are assigned global (actually subdomain) indices in a node major fashion: The dofs of the first node are 
    /// numbered, then the dofs of the second node, etc. Contrary to <see cref="NodeMajorDofOrderer"/>, the dofs of each node
    /// are assumed to be the same and supplied by the client. Based on that assumption, this class is much faster than
    /// <see cref="NodeMajorDofOrderer"/>. Constrained dofs are ignored.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class UniformDofOrderer: FreeDofOrdererBase
    {
        private readonly IReadOnlyList<DOFType> dofsPerNode;

        public UniformDofOrderer(IReadOnlyList<DOFType> dofsPerNode)
        {
            this.dofsPerNode = dofsPerNode;
            AreDofsOrdered = true;
        }

        public override void OrderDofs(ISubdomain subdomain)
        {
            FreeDofs = new DofTable();
            int dofCounter = 0;
            foreach (INode node in subdomain.Nodes)
            {
                bool isNodeConstrained = subdomain.Constraints.TryGetValue(node.ID, 
                    out Dictionary<DOFType, double> constraintsOfNode);
                foreach (DOFType dof in dofsPerNode)
                {
                    bool isDofConstrained = isNodeConstrained ? constraintsOfNode.ContainsKey(dof) : false;
                    if (!isDofConstrained) FreeDofs[node, dof] = dofCounter++;
                }
            }
            NumFreeDofs = dofCounter;
        }
    }
}
