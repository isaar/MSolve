using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcpg;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Projection
{
    internal class PreconditionerAsMatrixQ : IMatrixQ
    {
        private readonly IFetiPreconditioner preconditioner;

        internal PreconditionerAsMatrixQ(IFetiPreconditioner preconditioner)
        {
            this.preconditioner = preconditioner;
        }

        public Vector Multiply(Vector vector)
        {
            var result = Vector.CreateZero(vector.Length);
            preconditioner.SolveLinearSystem(vector, result);
            return result;
        }

        public Matrix Multiply(Matrix matrix)
        {
            var result = Matrix.CreateZero(matrix.NumRows, matrix.NumColumns);
            preconditioner.SolveLinearSystems(matrix, result);
            return result;
        }
    }
}
