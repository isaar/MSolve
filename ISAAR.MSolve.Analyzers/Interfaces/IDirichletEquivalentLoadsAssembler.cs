using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: _v2 delete the original one (IEquivalentLoadsAssembler) in FEM.Interfaces
namespace ISAAR.MSolve.Analyzers
{
    public interface IDirichletEquivalentLoadsAssembler
    {
        IVector GetEquivalentNodalLoads(ISubdomain_v2 subdomain, IVectorView solution, double constraintScalingFactor);
    }
}
