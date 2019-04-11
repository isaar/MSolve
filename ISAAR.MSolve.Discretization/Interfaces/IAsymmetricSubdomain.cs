using ISAAR.MSolve.Discretization.FreedomDegrees;

namespace ISAAR.MSolve.Discretization.Interfaces
{
	public interface IAsymmetricSubdomain:ISubdomain_v2
	{
		ISubdomainFreeDofOrdering DofRowOrdering { get; set; }
		ISubdomainFreeDofOrdering DofColOrdering { get; set; }
	}
}