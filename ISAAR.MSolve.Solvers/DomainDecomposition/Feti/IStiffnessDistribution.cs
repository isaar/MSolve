using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.Commons;

//TODO: This should be an enum class. There are only 2 possible cases.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    public interface IStiffnessDistribution
    {
        INodalLoadDistributor NodalLoadDistributor {get;}

        Dictionary<int, Matrix> CalcBoundaryPreconditioningSignedBooleanMatrices(DofSeparator dofSeparator, 
            LagrangeMultipliersEnumerator lagrangeEnumerator, Dictionary<int, Matrix> boundarySignedBooleanMatrices);

        //TODO: This only makes sense for HeterogeneousStiffnessDistribution. However, in the current design I must use the 
        //      IStiffnessDistribution in Feti1Solver. This constrained should be fixed though and this method removed from here.
        void StoreStiffnesses(Dictionary<int, IMatrixView> stiffnessMatrices,
            Table<INode, DOFType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses);
    }
}
