using System;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.Ordering.Reordering
{
    /// <summary>
    /// Does not apply any reordering. No object is mutated due to this class.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class NullReordering : IDofReorderingStrategy
    {
        public void ReorderDofs(ISubdomain_v2 subdomain, ISubdomainFreeDofOrdering originalOrdering)
        {
            // Do nothing
        }
    }
}
