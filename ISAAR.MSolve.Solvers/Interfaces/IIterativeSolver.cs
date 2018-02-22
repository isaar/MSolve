using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface IIterativeSolver : ISolver
    {
        int CurrentIteration { get; }
        void Initialize(IVector x, IVector residual, double detf);
        void Solve(int maxIterations, double tolerance);
    }
}
