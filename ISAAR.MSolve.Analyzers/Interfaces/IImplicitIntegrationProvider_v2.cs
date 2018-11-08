using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;
using System.Collections.Generic;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IImplicitIntegrationProvider_v2 : IAnalyzerProvider
    {
        void CalculateEffectiveMatrix(ILinearSystem_v2 subdomain, ImplicitIntegrationCoefficients coefficients);
        void ProcessRHS(ILinearSystem_v2 subdomain, ImplicitIntegrationCoefficients coefficients);
        IDictionary<int, Vector> GetAccelerationsOfTimeStep(int timeStep);
        IDictionary<int, Vector> GetVelocitiesOfTimeStep(int timeStep);
        void GetRHSFromHistoryLoad(int timeStep);
        IVector MassMatrixVectorProduct(ILinearSystem_v2 linearSystem, IVectorView lhsVector);
        IVector DampingMatrixVectorProduct(ILinearSystem_v2 linearSystem, IVectorView lhsVector);
    }
}
