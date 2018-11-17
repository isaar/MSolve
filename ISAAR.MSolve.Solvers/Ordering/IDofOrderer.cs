using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.Ordering
{
    public interface IDofOrderer
    {
        IDofOrdering OrderDofs(ISubdomain_v2 subdomain);
    }
}
