using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.LinearSystems
{
    /// <summary>
    /// Wrapper for a matrix class, so that it can be used by an iterative algorithm.
    /// </summary>
    public class MatrixTransformation: ILinearTransformation<Vector>
    {
        private readonly IMatrixView matrix;

        public MatrixTransformation(IMatrixView matrix)
        {
            this.matrix = matrix;
        }

        public Vector Multiply(Vector vector)
        {
            return matrix.MultiplyRight(vector, false);
        }
    }
}
