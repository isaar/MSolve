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
    public class NodeMajorDofOrderer: DofOrdererBase
    {
        protected override (int numGlobalFreeDofs, DofTable globalFreeDofs) OrderGlobalDofs(IStructuralModel_v2 model)
        {
            (int numGlobalFreeDofs, DofTable globalFreeDofs) = 
                SimpleDofOrderer.OrderFreeDofsOfElementSet(model.Elements, model.Constraints);
            globalFreeDofs.ReorderNodeMajor(model.Nodes);
            return (numGlobalFreeDofs, globalFreeDofs);
        }

        protected override (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) OrderSubdomainDofs(ISubdomain_v2 subdomain)
        {
            (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) =
                    SimpleDofOrderer.OrderFreeDofsOfElementSet(subdomain.Elements, subdomain.Constraints);
            subdomainFreeDofs.ReorderNodeMajor(subdomain.Nodes);
            return (numSubdomainFreeDofs, subdomainFreeDofs);
        }
    }
}
