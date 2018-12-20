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
        private readonly bool doOptimizationsIfSingleSubdomain = true;
        private readonly IDofOrderingStrategy orderingStrategy;
        private readonly IDofReorderingStrategy reorderingStrategy;

        public DofOrderer(IDofOrderingStrategy orderingStrategy, IDofReorderingStrategy reorderingStrategy,
            bool doOptimizationsIfSingleSubdomain = true, bool cacheElementToSubdomainDofMaps = true)
        {
            this.orderingStrategy = orderingStrategy;
            this.reorderingStrategy = reorderingStrategy;
            this.doOptimizationsIfSingleSubdomain = doOptimizationsIfSingleSubdomain;
            this.cacheElementToSubdomainDofMaps = cacheElementToSubdomainDofMaps;
        }

        public IReorderingAlgorithm Reordering { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        public IGlobalFreeDofOrdering OrderDofs(IStructuralModel_v2 model)
        {
            if (doOptimizationsIfSingleSubdomain && (model.Subdomains.Count == 1))
            {
                ISubdomain_v2 subdomain = model.Subdomains.First();

                // Order subdomain dofs
                (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) = orderingStrategy.OrderSubdomainDofs(subdomain);
                ISubdomainFreeDofOrdering subdomainOrdering;
                if (cacheElementToSubdomainDofMaps) subdomainOrdering = new SubdomainFreeDofOrderingCaching(
                    numSubdomainFreeDofs, subdomainFreeDofs);
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
                    (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) = orderingStrategy.OrderSubdomainDofs(subdomain);
                    ISubdomainFreeDofOrdering subdomainOrdering;
                    if (cacheElementToSubdomainDofMaps) subdomainOrdering = new SubdomainFreeDofOrderingCaching(
                        numSubdomainFreeDofs, subdomainFreeDofs);
                    else subdomainOrdering = new SubdomainFreeDofOrderingGeneral(numSubdomainFreeDofs, subdomainFreeDofs);
                    subdomainOrderings.Add(subdomain, subdomainOrdering);

                    // Reorder subdomain dofs
                    reorderingStrategy.ReorderDofs(subdomain, subdomainOrdering);
                }

                // Order global dofs
                (int numGlobalFreeDofs, DofTable globalFreeDofs) = orderingStrategy.OrderGlobalDofs(model);
                return new GlobalFreeDofOrderingGeneral(numGlobalFreeDofs, globalFreeDofs, subdomainOrderings);
            }
        }
    }
}
