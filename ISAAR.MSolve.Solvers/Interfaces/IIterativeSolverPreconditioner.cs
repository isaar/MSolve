using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface IIterativeSolverPreconditioner
    {
        void Precondition(IVectorOLD vIn, IVectorOLD vOut);
    }
}
