using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;
using System.Collections.Generic;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IImplicitIntegrationProvider_v2 : IAnalyzerProvider
    {
        void CalculateEffectiveMatrix(ILinearSystem_v2 subdomain, ImplicitIntegrationCoefficients coefficients);
        void ProcessRhs(ILinearSystem_v2 subdomain, ImplicitIntegrationCoefficients coefficients);
        IDictionary<int, IVector> GetAccelerationsOfTimeStep(int timeStep);
        IDictionary<int, IVector> GetVelocitiesOfTimeStep(int timeStep);
        void GetRhsFromHistoryLoad(int timeStep);

        //TODO: what about thermal? There is no mass matrix there. Either define these as 1st order matrix coeff, 2nd order ...
        //      or let the provider call them instead of the analyzer.
        IVector MassMatrixVectorProduct(ILinearSystem_v2 linearSystem, IVectorView lhsVector); 
        IVector DampingMatrixVectorProduct(ILinearSystem_v2 linearSystem, IVectorView lhsVector);
    }
}
