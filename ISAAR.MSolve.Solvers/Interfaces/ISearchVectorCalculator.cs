using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface ISearchVectorCalculator
    {
        bool InitializeStartingVectorFromSearchVectors(IVectorOLD x, IVectorOLD b);
        void CalculateSearchVector(IIterativeSolver solver);
        double CalculateGradient(IIterativeSolver solver);
        void ClearSearchVectors(int vectorsToKeepFromTop);
    }
}
