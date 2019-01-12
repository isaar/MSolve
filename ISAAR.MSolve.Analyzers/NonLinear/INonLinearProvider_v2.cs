using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Analyzers.NonLinear
{
    public interface INonLinearProvider_v2 : IAnalyzerProvider_v2
    {
        double CalculateRhsNorm(IVectorView rhs);

        //TODO: Very generic name. There is also a similar method in IImplictIntegrationProvider.
        void ProcessInternalRhs(ISubdomain_v2 subdomain, IVectorView solution, IVector rhs);
    }
}
