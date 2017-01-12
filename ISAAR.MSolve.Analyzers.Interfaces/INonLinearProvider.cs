using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface INonLinearProvider : IAnalyzerProvider
    {
        double RHSNorm(double[] rhs);
        void ProcessInternalRHS(IMatrix2D matrix, double[] rhs, double[] solution);
    }
}
