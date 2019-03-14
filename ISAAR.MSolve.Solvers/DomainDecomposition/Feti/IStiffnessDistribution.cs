using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: This should be an enum class. There are only 2 possible cases.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    public interface IStiffnessDistribution
    {
        INodalLoadDistributor NodalLoadDistributor {get;}

        //TODO: Pass in an object that stores dof data, including the multiplicities
        Matrix CalcBoundaryPreconditioningSignedBooleanMatrix(Matrix boundarySignedBooleanMatrix,
            int[] boundaryDofsMultiplicity);
    }
}
