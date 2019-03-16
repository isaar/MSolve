using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: During preconditioning (including operation with the Q matrix), each subdomain performs multiplications with 
//      Bpb^T*Sbb*Bpb * x. This matrix could be explicitly formed and cached
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    public interface IFetiPreconditioner
    {
        //TODO: Would it be better to just return the results instead of overwriting them? It would be easier to use and no one would have to do Clear().
        void SolveLinearSystem(Vector rhs, Vector lhs);
        void SolveLinearSystems(Matrix rhs, Matrix lhs);
    }

    public interface IFetiPreconditionerFactory
    {
        IFetiPreconditioner CreatePreconditioner(IStiffnessDistribution stiffnessDistribution, DofSeparator dofSeparator, 
            LagrangeMultipliersEnumerator lagrangeEnumerator, Dictionary<int, IMatrixView> stiffnessMatrices);
    }
}
