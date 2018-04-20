using ISAAR.MSolve.Solvers.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IStaticProvider : IAnalyzerProvider
    {
        void CalculateMatrix(ILinearSystem subdomain);
        //void CalculateMatrices();
    }
}
