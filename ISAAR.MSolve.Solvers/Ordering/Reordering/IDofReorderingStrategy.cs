using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.Ordering.Reordering
{
    /// <summary>
    /// Reorders the freedom degrees of a subdomain.
    /// Author: Serafeim Bakalakos
    /// </summary>
    public interface IDofReorderingStrategy
    {
        void ReorderDofs(ISubdomain_v2 subdomain, ISubdomainFreeDofOrdering originalOrdering);
    }
}
