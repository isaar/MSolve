using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface ISearchVectorCalculator
    {
        bool InitializeStartingVectorFromSearchVectors(IVector x, IVector b);
        void CalculateSearchVector(IIterativeSolver solver);
        double CalculateGradient(IIterativeSolver solver);
        void ClearSearchVectors(int vectorsToKeepFromTop);
    }
}
