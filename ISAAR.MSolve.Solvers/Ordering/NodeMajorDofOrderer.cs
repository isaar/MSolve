using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// Free dofs are assigned global (actually subdomain) indices in a node major fashion: The dofs of the first node are 
    /// numbered, then the dofs of the second node, etc. Constrained dofs are ignored.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class NodeMajorDofOrderer: FreeDofOrdererBase
    {
        public override void OrderDofs(ISubdomain subdomain)
        {
            (NumFreeDofs, FreeDofs) = SimpleDofOrderer.OrderFreeDofsAtFirstOccurence(subdomain);
            FreeDofs.ReorderNodeMajor(subdomain.Nodes);
            AreDofsOrdered = true;
        }
    }
}
