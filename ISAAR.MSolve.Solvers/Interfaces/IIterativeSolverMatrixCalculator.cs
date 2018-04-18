using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface IIterativeSolverMatrixCalculator
    {
        void MultiplyWithMatrix(IVector vIn, IVector vOut);
    }
}
