using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// Free dofs are assigned global (actually subdomain) indices in a node major fashion: The dofs of the first node are 
    /// numbered, then the dofs of the second node, etc. Constrained dofs are ignored.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class NodeMajorDofOrderer: IDofOrderer
    {
        private readonly SimpleDofOrderer embeddedOrderer = new SimpleDofOrderer();

        public IDofOrdering OrderDofs(ISubdomain subdomain)
        {
            IDofOrdering ordering = embeddedOrderer.OrderDofs(subdomain);
            ordering.FreeDofs.ReorderNodeMajor(subdomain.Nodes);
            return ordering;
        }
    }
}
