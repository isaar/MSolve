using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: This should be called second order provider. The matrices, coefficients, etc. should be named 0th-order, 1st-order, 
//      2nd-order. 
//TODO: Implicit/explicit time integration logic should be defined by the analyzer and implemented by the provider, in order to 
//      reuse the analyzer for problems that have a slightly different differential equation (e.g. coupled problems).
//TODO: Perhaps the providers should not hold references to the linear systems. Instead they would return vectors/matrices to 
//      the analyzers (or the vectors/matrices would be passed in and overwritten).
//TODO: Rename the Get~ methods to Calculate or something similar.
namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IImplicitIntegrationProvider : IAnalyzerProvider
    {
        //TODO: This should not exist at all. The provider should return the 0th order (stiffness), 1st order (damping) and 2nd
        //      order matrices (or some matrix representations that can be combined between them and multiplied with vectors).
        IMatrixView LinearCombinationOfMatricesIntoStiffness(ImplicitIntegrationCoefficients coefficients, 
            ISubdomain subdomain);

        //TODO: Way too generic name. Probably needs refactoring as well.
        void ProcessRhs(ImplicitIntegrationCoefficients coefficients, ISubdomain subdomain, IVector rhs);

        //TODO: I think the analyzer is responsible for these. E.g. Newmark has the formulas with beta and gamma, Euler has the
        //      central differences formulas, etc.
        IDictionary<int, IVector> GetAccelerationsOfTimeStep(int timeStep);
        IDictionary<int, IVector> GetVelocitiesOfTimeStep(int timeStep);

        IDictionary<int, IVector> GetRhsFromHistoryLoad(int timeStep);

        //TODO: what about thermal? There is no mass matrix there. Define these as 1st order matrix coeff, 2nd order ...
        IVector MassMatrixVectorProduct(ISubdomain subdomain, IVectorView vector); 
        IVector DampingMatrixVectorProduct(ISubdomain subdomain, IVectorView vector);
    }
}
