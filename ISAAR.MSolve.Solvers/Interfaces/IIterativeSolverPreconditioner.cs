using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface IIterativeSolverPreconditioner
    {
        void Precondition(IVector vIn, IVector vOut);
    }
}
