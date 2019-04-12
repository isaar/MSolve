using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface INonLinearProvider : IAnalyzerProvider
    {
        double CalculateRhsNorm(IVectorView rhs);

        //TODO: Very generic name. There is also a similar method in IImplictIntegrationProvider.
        void ProcessInternalRhs(ISubdomain subdomain, IVectorView solution, IVector rhs);
    }
}
