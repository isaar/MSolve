using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Solvers.PCG
{
    public class JacobiPreconditionerBuilder : IPreconditionerBuilder
    {
        public IPreconditioner BuildPreconditioner(IMatrixView matrix)
        {
            return new JacobiPreconditioner(matrix.GetDiagonalAsArray());
        }
    }
}
