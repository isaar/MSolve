using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IImplicitIntegrationProvider : IAnalyzerProvider
    {
        void CalculateEffectiveMatrix(IMatrix2D matrix, IImplicitIntegrationCoefficients coefficients);
        void ProcessRHS(IMatrix2D matrix, IImplicitIntegrationCoefficients coefficients);
        void GetRHSFromHistoryLoad(int timeStep);
        void MassMatrixVectorProduct(IMatrix2D matrix, IVector vIn, double[] vOut);
        void DampingMatrixVectorProduct(IMatrix2D matrix, IVector vIn, double[] vOut);
    }
}
