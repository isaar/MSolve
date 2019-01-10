using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: Not sure if this is needed.
namespace ISAAR.MSolve.Solvers.Ordering
{
    public interface IDofOrderer
    {
        IGlobalFreeDofOrdering OrderDofs(IStructuralModel_v2 model);
    }
}
