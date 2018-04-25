using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using System.Collections.Generic;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IImplicitIntegrationProvider : IAnalyzerProvider
    {
        void CalculateEffectiveMatrix(ILinearSystem subdomain, ImplicitIntegrationCoefficients coefficients);
        void ProcessRHS(ILinearSystem subdomain, ImplicitIntegrationCoefficients coefficients);
        IDictionary<int, double[]> GetAccelerationsOfTimeStep(int timeStep);
        IDictionary<int, double[]> GetVelocitiesOfTimeStep(int timeStep);
        void GetRHSFromHistoryLoad(int timeStep);
        void MassMatrixVectorProduct(ILinearSystem subdomain, IVector vIn, double[] vOut);
        void DampingMatrixVectorProduct(ILinearSystem subdomain, IVector vIn, double[] vOut);
    }
}
