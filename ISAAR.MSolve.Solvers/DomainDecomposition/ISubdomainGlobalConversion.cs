using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: the operations described here are expressed as matrix-vector multiplications with the projection & scaling matrices 
//      called Lpb in FETI bibliography. Should I also express them similarly? In that case dedicated classes are needed for 
//      the Lpb matrices, but the code for the conversions is the same for homogeneous/heterogeneous. Only the creation of the
//      Lpb matrices depends on the stiffness distribution, which is expected and needed elsewhere too.
//TODO: These methods overlap with IGlobalDofOrdering. Perhaps they should be moved there, while the implementing classes remain 
//      here as strategies.
//TODO: Should the vectors be Vector instead of IVectorView?
namespace ISAAR.MSolve.Solvers.DomainDecomposition
{
    public interface ISubdomainGlobalConversion
    {
        double CalculateGlobalForcesNorm(Dictionary<int, IVectorView> subdomainForces);

        Dictionary<int, SparseVector> DistributeNodalLoads(IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems,
            Table<INode, DOFType, double> globalNodalLoads);

        //TODO: arguments other than the displacements should be injected into the constructor.
        Vector GatherGlobalDisplacements(Dictionary<int, IVectorView> subdomainDisplacements);

        Vector GatherGlobalForces(Dictionary<int, IVectorView> subdomainForces);
    }
}
