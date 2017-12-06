using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IImplicitIntegrationProvider : IAnalyzerProvider
    {
        // TODO: ID fields should be removed and the whole domain decomposition thing must be decoupled from the analyzers
        IMatrix2D CalculateEffectiveMatrix(int id, IImplicitIntegrationCoefficients coefficients);
        void ProcessRHS(int id, IImplicitIntegrationCoefficients coefficients);
        void GetRHSFromHistoryLoad(int timeStep);
        void MassMatrixVectorProduct(int id, IVectorOLD vIn, double[] vOut);
        void DampingMatrixVectorProduct(int id, IVectorOLD vIn, double[] vOut);
    }
}
