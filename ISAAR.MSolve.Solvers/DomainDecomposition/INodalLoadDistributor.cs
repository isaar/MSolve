using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: This should be an enum class. There are only 2 possible cases.
//TODO: Perhaps it should be incorporated into IStiffnessDistribution. Otherwise calculating rhs norms should be done here.
namespace ISAAR.MSolve.Solvers.DomainDecomposition
{
    public interface INodalLoadDistributor
    {
        Dictionary<int, SparseVector> DistributeNodalLoads(IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems,
            Table<INode, DOFType, double> globalNodalLoads);
    }
}
