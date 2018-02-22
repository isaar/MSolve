using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IImplicitIntegrationProvider : IAnalyzerProvider
    {
        void CalculateEffectiveMatrix(ILinearSystem subdomain, ImplicitIntegrationCoefficients coefficients);
        void ProcessRHS(ILinearSystem subdomain, ImplicitIntegrationCoefficients coefficients);
        void GetRHSFromHistoryLoad(int timeStep);
        void MassMatrixVectorProduct(ILinearSystem subdomain, IVector vIn, double[] vOut);
        void DampingMatrixVectorProduct(ILinearSystem subdomain, IVector vIn, double[] vOut);
    }
}
