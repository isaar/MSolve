using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using System.Collections.Generic;

//TODO: This should be called second order provider. Implicit/explicit time integration logic should be handled by the 
//      appropriate analyzer, not the providers.
//TODO: Perhaps the providers should not hold references to the linear systems. Instead they would return vectors/matrices to 
//      the analyzers (or the vectors/matrices would be passed in and overwritten).
namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IImplicitIntegrationProvider_v2 : IAnalyzerProvider
    {
        //TODO: This should not exist at all. The provider should return the 0th order (stiffness), 1st order (damping) and 2nd
        //      order matrices (or some matrix representations that can be combined between them and multiplied with vectors).
        void CalculateEffectiveMatrix(ILinearSystem_v2 subdomain, ImplicitIntegrationCoefficients coefficients);

        //TODO: Way too generic name. Probably needs refactoring as well.
        void ProcessRhs(ILinearSystem_v2 subdomain, ImplicitIntegrationCoefficients coefficients);

        //TODO: I think the analyzer is responsible for these. E.g. Newmark has the formulas with beta and gamma, Euler has the
        //      central differences formulas, etc.
        IDictionary<int, IVector> GetAccelerationsOfTimeStep(int timeStep);
        IDictionary<int, IVector> GetVelocitiesOfTimeStep(int timeStep);

        void GetRhsFromHistoryLoad(int timeStep);

        //TODO: what about thermal? There is no mass matrix there. Define these as 1st order matrix coeff, 2nd order ...
        IVector MassMatrixVectorProduct(ILinearSystem_v2 linearSystem, IVectorView lhsVector); 
        IVector DampingMatrixVectorProduct(ILinearSystem_v2 linearSystem, IVectorView lhsVector);
    }
}
