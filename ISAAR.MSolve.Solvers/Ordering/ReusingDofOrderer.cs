using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: This should be named DofOrderer and the old one should be removed. I kept the old version for backwards compatibility 
//      while developing & testing this
namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// Orders the unconstrained freedom degrees of each subdomain and the shole model. Also applies any reordering and other 
    /// optimizations.
    /// </summary>
    public class ReusingDofOrderer : IDofOrderer
    {
        //TODO: this should also be a strategy, so that I could have caching with fallbacks, in case of insufficient memor.
        private readonly bool cacheElementToSubdomainDofMaps = true;
        private readonly bool doOptimizationsIfSingleSubdomain = true; // No idea why someone would want this to be false.
        private readonly ConstrainedDofOrderingStrategy constrainedOrderingStrategy;
        private readonly IFreeDofOrderingStrategy freeOrderingStrategy;
        private readonly IDofReorderingStrategy reorderingStrategy;

        public ReusingDofOrderer(IFreeDofOrderingStrategy freeOrderingStrategy, IDofReorderingStrategy reorderingStrategy,
            bool doOptimizationsIfSingleSubdomain = true, bool cacheElementToSubdomainDofMaps = true)
        {
            this.constrainedOrderingStrategy = new ConstrainedDofOrderingStrategy();
            this.freeOrderingStrategy = freeOrderingStrategy;
            this.reorderingStrategy = reorderingStrategy;
            this.doOptimizationsIfSingleSubdomain = doOptimizationsIfSingleSubdomain;
            this.cacheElementToSubdomainDofMaps = cacheElementToSubdomainDofMaps;
        }

        public ISubdomainConstrainedDofOrdering OrderConstrainedDofs(ISubdomain subdomain)
        {
            (int numConstrainedDofs, DofTable constrainedDofs) =
                constrainedOrderingStrategy.OrderSubdomainDofs(subdomain);
            if (cacheElementToSubdomainDofMaps)
            {
                return new SubdomainConstrainedDofOrderingCaching(numConstrainedDofs, constrainedDofs);
            }
            else return new SubdomainConstrainedDofOrderingGeneral(numConstrainedDofs, constrainedDofs);
        }

        public IGlobalFreeDofOrdering OrderFreeDofs(IStructuralModel model)
        {
            if (doOptimizationsIfSingleSubdomain && (model.Subdomains.Count == 1))
            {
                ISubdomain subdomain = model.Subdomains.First();
                ISubdomainFreeDofOrdering subdomainOrdering;
                if (!subdomain.ConnectivityModified) subdomainOrdering = subdomain.FreeDofOrdering;
                else
                {
                    // Order subdomain dofs
                    Debug.WriteLine($"{this.GetType().Name}: Ordering free dofs of subdomain {subdomain.ID}");
                    (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) = freeOrderingStrategy.OrderSubdomainDofs(subdomain);
                    if (cacheElementToSubdomainDofMaps)
                    {
                        subdomainOrdering = new SubdomainFreeDofOrderingCaching(numSubdomainFreeDofs, subdomainFreeDofs);
                    }
                    else subdomainOrdering = new SubdomainFreeDofOrderingGeneral(numSubdomainFreeDofs, subdomainFreeDofs);

                    // Reorder subdomain dofs
                    reorderingStrategy.ReorderDofs(subdomain, subdomainOrdering);
                }
                
                // Order global dofs
                return new GlobalFreeDofOrderingSingle(subdomain, subdomainOrdering);
            }
            else
            {
                // Order subdomain dofs
                var subdomainOrderings = new Dictionary<ISubdomain, ISubdomainFreeDofOrdering>(model.Subdomains.Count);
                foreach (ISubdomain subdomain in model.Subdomains)
                {
                    ISubdomainFreeDofOrdering subdomainOrdering;
                    if (!subdomain.ConnectivityModified) subdomainOrdering = subdomain.FreeDofOrdering;
                    else 
                    {
                        Debug.WriteLine($"{this.GetType().Name}: Ordering free dofs of subdomain {subdomain.ID}");
                        (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) = freeOrderingStrategy.OrderSubdomainDofs(subdomain);
                        if (cacheElementToSubdomainDofMaps)
                        {
                            subdomainOrdering = new SubdomainFreeDofOrderingCaching(numSubdomainFreeDofs, subdomainFreeDofs);

                        }
                        else subdomainOrdering = new SubdomainFreeDofOrderingGeneral(numSubdomainFreeDofs, subdomainFreeDofs);

                        // Reorder subdomain dofs
                        reorderingStrategy.ReorderDofs(subdomain, subdomainOrdering);
                    }
                    subdomainOrderings.Add(subdomain, subdomainOrdering);
                }

                // Order global dofs
                (int numGlobalFreeDofs, DofTable globalFreeDofs) = freeOrderingStrategy.OrderGlobalDofs(model);
                return new GlobalFreeDofOrderingGeneral(numGlobalFreeDofs, globalFreeDofs, subdomainOrderings);
            }
        }
    }
}
