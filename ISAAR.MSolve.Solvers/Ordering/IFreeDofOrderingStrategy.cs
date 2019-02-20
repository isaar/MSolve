using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// Determines how the unconstrained freedom degrees of the physical model will be ordered.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IFreeDofOrderingStrategy
    {
        /// <summary>
        /// Orders the unconstrained freedom degrees of the whole model.
        /// </summary>
        /// <param name="model">The physical model that is analyzed.</param>
        (int numGlobalFreeDofs, DofTable globalFreeDofs) OrderGlobalDofs(IStructuralModel_v2 model);

        /// <summary>
        /// Orders the unconstrained freedom degrees of one of the model's subdomains.
        /// </summary>
        /// <param name="subdomain">A subdomain of the whole model.</param>
        (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) OrderSubdomainDofs(ISubdomain_v2 subdomain);
    }
}
