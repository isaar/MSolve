using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.FETI
{
    internal interface IFetiPreconditioner
    {
        void CreatePreconditioner(Dictionary<int, IMatrixView> stiffnessMatrices);
        void SolveLinearSystem(Vector rhs, Vector lhs);
    }
}
