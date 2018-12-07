using ISAAR.MSolve.Solvers.Commons;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IStaticProvider_v2 : IAnalyzerProvider
    {
        //TODO: this should calculate and return the matrix and the analyzer should set it. Perhaps that convention should hold 
        //      for all matrices, rhs vectors.
        void CalculateMatrix(ILinearSystem_v2 subdomain);
    }
}
