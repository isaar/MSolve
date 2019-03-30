using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// Orders the unconstrained freedom degrees of each subdomain and the shole model. Also applies any reordering and other 
    /// optimizations.
    /// </summary>
    public class DofOrderer : IDofOrderer
    {
        //TODO: this should also be a strategy, so that I could have caching with fallbacks, in case of insufficient memor.
        private readonly bool cacheElementToSubdomainDofMaps = true; 
        private readonly bool doOptimizationsIfSingleSubdomain = true; // No idea why someone would want this to be false.
        private readonly ConstrainedDofOrderingStrategy constrainedOrderingStrategy;
        private readonly IFreeDofOrderingStrategy freeOrderingStrategy;
        private readonly IDofReorderingStrategy reorderingStrategy;

        public DofOrderer(IFreeDofOrderingStrategy freeOrderingStrategy, IDofReorderingStrategy reorderingStrategy,
            bool doOptimizationsIfSingleSubdomain = true, bool cacheElementToSubdomainDofMaps = true)
        {
            this.constrainedOrderingStrategy = new ConstrainedDofOrderingStrategy();
            this.freeOrderingStrategy = freeOrderingStrategy;
            this.reorderingStrategy = reorderingStrategy;
            this.doOptimizationsIfSingleSubdomain = doOptimizationsIfSingleSubdomain;
            this.cacheElementToSubdomainDofMaps = cacheElementToSubdomainDofMaps;
        }

        public ISubdomainConstrainedDofOrdering OrderConstrainedDofs(ISubdomain_v2 subdomain)
        {
            (int numConstrainedDofs, DofTable constrainedDofs) = 
                constrainedOrderingStrategy.OrderSubdomainDofs(subdomain);
            if (cacheElementToSubdomainDofMaps)
            {
                return new SubdomainConstrainedDofOrderingCaching(numConstrainedDofs, constrainedDofs);
            }
            else return new SubdomainConstrainedDofOrderingGeneral(numConstrainedDofs, constrainedDofs);
        }

        public IGlobalFreeDofOrdering OrderFreeDofs(IStructuralModel_v2 model)
        {
            if (doOptimizationsIfSingleSubdomain && (model.Subdomains.Count == 1))
            {
                ISubdomain_v2 subdomain = model.Subdomains.First();

                // Order subdomain dofs
                (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) = freeOrderingStrategy.OrderSubdomainDofs(subdomain);
                ISubdomainFreeDofOrdering subdomainOrdering;
                if (cacheElementToSubdomainDofMaps)
                {
                    subdomainOrdering = new SubdomainFreeDofOrderingCaching(numSubdomainFreeDofs, subdomainFreeDofs);
                }
                else subdomainOrdering = new SubdomainFreeDofOrderingGeneral(numSubdomainFreeDofs, subdomainFreeDofs);

                // Reorder subdomain dofs
                reorderingStrategy.ReorderDofs(subdomain, subdomainOrdering);

                // Order global dofs
                return new GlobalFreeDofOrderingSingle(subdomain, subdomainOrdering);
            }
            else
            {
                // Order subdomain dofs
                var subdomainOrderings = new Dictionary<ISubdomain_v2, ISubdomainFreeDofOrdering>(model.Subdomains.Count);
                foreach (ISubdomain_v2 subdomain in model.Subdomains)
                {
                    (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) = freeOrderingStrategy.OrderSubdomainDofs(subdomain);
                    ISubdomainFreeDofOrdering subdomainOrdering;
                    if (cacheElementToSubdomainDofMaps) subdomainOrdering = new SubdomainFreeDofOrderingCaching(
                        numSubdomainFreeDofs, subdomainFreeDofs);
                    else subdomainOrdering = new SubdomainFreeDofOrderingGeneral(numSubdomainFreeDofs, subdomainFreeDofs);
                    subdomainOrderings.Add(subdomain, subdomainOrdering);

                    // Reorder subdomain dofs
                    reorderingStrategy.ReorderDofs(subdomain, subdomainOrdering);
                }

                // Order global dofs
                (int numGlobalFreeDofs, DofTable globalFreeDofs) = freeOrderingStrategy.OrderGlobalDofs(model);
                return new GlobalFreeDofOrderingGeneral(numGlobalFreeDofs, globalFreeDofs, subdomainOrderings);
            }
        }
    }
}
