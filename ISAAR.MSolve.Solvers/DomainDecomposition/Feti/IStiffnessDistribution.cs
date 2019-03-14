using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;

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
