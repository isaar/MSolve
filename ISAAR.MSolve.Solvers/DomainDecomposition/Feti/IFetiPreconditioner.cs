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
        IFetiPreconditioner CreatePreconditioner(Dictionary<int, int[]> boundaryDofs, Dictionary<int, int[]> internalDofs,
            ContinuityEquationsCalculator continuityEquations, Dictionary<int, IMatrixView> stiffnessMatrices);
    }
}
