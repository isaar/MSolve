using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.Ordering
{
    public interface IDofOrderer
    {
        IDofOrdering OrderDofs(IStructuralModel_v2 model, ISubdomain_v2 subdomain);
    }
}
