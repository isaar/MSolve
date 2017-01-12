using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface INonLinearParentAnalyzer : IAnalyzer
    {
        double[] GetOtherRHSComponents(IMatrix2D matrix, IVector currentSolution);
    }
}
