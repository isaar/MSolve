using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IStaticProvider_v2 : IAnalyzerProvider_v2
    {
        IMatrixView CalculateMatrix(ISubdomain_v2 subdomain);
    }
}
