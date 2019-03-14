using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    public interface IFetiPreconditioner
    {
        void SolveLinearSystem(Vector rhs, Vector lhs);
    }

    public interface IFetiPreconditionerFactory
    {
        //TODO: dof data should be calculated and managed by a single class
        IFetiPreconditioner CreatePreconditioner(IStiffnessDistribution stiffnessDistribution,
            Dictionary<int, int[]> boundaryDofs, Dictionary<int, int[]> boundaryDofsMultiplicity,
            Dictionary<int, int[]> internalDofs, ContinuityEquationsCalculator continuityEquations, 
            Dictionary<int, IMatrixView> stiffnessMatrices);
    }
}
