using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface IIterativeSolver : ISolver
    {
        int CurrentIteration { get; }
        void Initialize(IVectorOLD x, IVectorOLD residual, double detf);
        void Solve(int maxIterations, double tolerance);
    }
}
