using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.Ordering.Reordering
{
    /// <summary>
    /// Reorders the unconstrained freedom degrees by using multiple other reordering strategies. For now, only consecutive 
    /// application of each reordering strategy is possible.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CompositeReordering : IDofReorderingStrategy
    {
        private readonly IReadOnlyList<IDofReorderingStrategy> reorderingStrategies;

        public CompositeReordering(params IDofReorderingStrategy[] reorderingStrategies)
        {
            this.reorderingStrategies = reorderingStrategies;
        }

        public void ReorderDofs(ISubdomain_v2 subdomain, ISubdomainFreeDofOrdering originalOrdering)
        {
            foreach (var reordering in reorderingStrategies) reordering.ReorderDofs(subdomain, originalOrdering);
        }
    }
}
