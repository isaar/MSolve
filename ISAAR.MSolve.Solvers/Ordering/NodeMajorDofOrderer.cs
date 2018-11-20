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

        public IGlobalFreeDofOrdering OrderDofs(IStructuralModel_v2 model)
        {
            //TODO: move this to the end
            (int numGlobalFreeDofs, DofTable globalFreeDofs) =
                   SimpleDofOrderer.OrderFreeDofsOfElementSet(model.Elements, model.Constraints);
            globalFreeDofs.ReorderNodeMajor(model.Nodes);

            // Order subdomain dofs
            var subdomainOrderings = new Dictionary<ISubdomain_v2, ISubdomainFreeDofOrdering>(model.Subdomains.Count);
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) =
                    SimpleDofOrderer.OrderFreeDofsOfElementSet(subdomain.Elements, subdomain.Constraints);
                subdomainFreeDofs.ReorderNodeMajor(subdomain.Nodes);
                ISubdomainFreeDofOrdering subdomainOrdering =
                    new SubdomainFreeDofOrderingGeneral(numSubdomainFreeDofs, subdomainFreeDofs, globalFreeDofs);
                subdomainOrderings.Add(subdomain, subdomainOrdering);
            }

            // Order global dofs
            return new GlobalFreeDofOrderingGeneral(numGlobalFreeDofs, globalFreeDofs, subdomainOrderings);
        }
    }
}
