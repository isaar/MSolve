using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

//TODO: This should be an enum class. There are only 2 possible cases.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.StiffnessDistribution
{
    public interface IFeti1StiffnessDistribution
    {
        //TODO: Ideally, ISubdomainGlobalConversion would be an independent class that only needs Lpb to operate and can be 
        //      reused in other solvers. IStiffnessDistribution would then be used only to create Lpb.
        ISubdomainGlobalConversion SubdomainGlobalConversion { get;}

        Dictionary<int, Matrix> CalcBoundaryPreconditioningSignedBooleanMatrices(
            LagrangeMultipliersEnumerator lagrangeEnumerator, Dictionary<int, Matrix> boundarySignedBooleanMatrices);
    }
}
