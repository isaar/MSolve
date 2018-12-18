using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Reordering;

namespace ISAAR.MSolve.Solvers.Ordering
{
    public interface IDofOrderer
    {
        IReorderingAlgorithm Reordering { get; set; }
        IGlobalFreeDofOrdering OrderDofs(IStructuralModel_v2 model);
    }
}
