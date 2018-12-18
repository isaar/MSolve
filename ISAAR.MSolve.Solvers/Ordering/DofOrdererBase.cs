using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Reordering;

//TODO: Is the injection of reordering done in the best way?
namespace ISAAR.MSolve.Solvers.Ordering
{
    public abstract class DofOrdererBase: IDofOrderer
    {
        public bool CacheElementToSubdomainDofMaps { get; set; } = true;
        public bool DoOptimizationsIfSingleSubdomain { get; set; } = true;

        public IReorderingAlgorithm Reordering { get; set; } = null; 

        public IGlobalFreeDofOrdering OrderDofs(IStructuralModel_v2 model)
        {
            if (DoOptimizationsIfSingleSubdomain && (model.Subdomains.Count == 1))
            {
                ISubdomain_v2 subdomain = model.Subdomains.First();

                // Order subdomain dofs
                (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) = OrderSubdomainDofs(subdomain);
                ISubdomainFreeDofOrdering subdomainOrdering;
                if (CacheElementToSubdomainDofMaps) subdomainOrdering = new SubdomainFreeDofOrderingCaching(
                    numSubdomainFreeDofs, subdomainFreeDofs);
                else subdomainOrdering = new SubdomainFreeDofOrderingGeneral(numSubdomainFreeDofs, subdomainFreeDofs);

                // Reorder subdomain dofs
                ReorderDofs(subdomain, subdomainOrdering);

                // Order global dofs
                return new GlobalFreeDofOrderingSingle(subdomain, subdomainOrdering);
            }
            else
            {
                // Order subdomain dofs
                var subdomainOrderings = new Dictionary<ISubdomain_v2, ISubdomainFreeDofOrdering>(model.Subdomains.Count);
                foreach (ISubdomain_v2 subdomain in model.Subdomains)
                {
                    (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) = OrderSubdomainDofs(subdomain);
                    ISubdomainFreeDofOrdering subdomainOrdering;
                    if (CacheElementToSubdomainDofMaps) subdomainOrdering = new SubdomainFreeDofOrderingCaching(
                        numSubdomainFreeDofs, subdomainFreeDofs);
                    else subdomainOrdering = new SubdomainFreeDofOrderingGeneral(numSubdomainFreeDofs, subdomainFreeDofs);
                    subdomainOrderings.Add(subdomain, subdomainOrdering);

                    // Reorder subdomain dofs
                    ReorderDofs(subdomain, subdomainOrdering);
                }

                // Order global dofs
                (int numGlobalFreeDofs, DofTable globalFreeDofs) = OrderGlobalDofs(model);
                return new GlobalFreeDofOrderingGeneral(numGlobalFreeDofs, globalFreeDofs, subdomainOrderings);
            }
        }

        protected void ReorderDofs(ISubdomain_v2 subdomain, ISubdomainFreeDofOrdering originalOrdering)
        {
            if (Reordering == null) return;
            var pattern = SparsityPatternSymmetric.CreateEmpty(originalOrdering.NumFreeDofs);
            foreach (var element in subdomain.Elements)
            {
                (int[] elementDofIndices, int[] subdomainDofIndices) = originalOrdering.MapFreeDofsElementToSubdomain(element);

                //TODO: ISubdomainFreeDofOrdering could perhaps return whether the subdomainDofIndices are sorted or not.
                pattern.ConnectIndices(subdomainDofIndices, false);
            }
            (int[] permutation, bool oldToNew) = Reordering.FindPermutation(pattern);
            originalOrdering.FreeDofs.Reorder(permutation, oldToNew);
        }

        protected abstract (int numGlobalFreeDofs, DofTable globalFreeDofs) OrderGlobalDofs(IStructuralModel_v2 model);
        protected abstract (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) OrderSubdomainDofs(ISubdomain_v2 subdomain);
    }
}
