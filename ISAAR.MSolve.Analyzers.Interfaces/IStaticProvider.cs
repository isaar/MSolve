using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IStaticProvider : IAnalyzerProvider
    {
        void CalculateMatrix(IMatrix2D subdomain);
    }
}
