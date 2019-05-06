using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.Assemblers.Collocation
{
	public interface IAsymmetricDofOrderingStrategy
	{
		(int numGlobalFreeDofs, DofTable globalFreeDofs) OrderGlobalDofs(IStructuralAsymmetricModel model);

		(int numSubdomainFreeDofs, DofTable subdomainFreeDofs) OrderSubdomainDofs(IAsymmetricSubdomain subdomain);
	}
}