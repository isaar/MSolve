using System;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// A diagonal matrix with arbitrary order and half its entries being negative. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class DiagonalIndefinite
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="evenOrder">The number of rows/columns. It must be even.</param>
        /// <returns></returns>
        internal static (Matrix A, Vector b, Vector x, IPreconditioner M) BuildIndefiniteSystem(int evenOrder)
        {
            // Example taken from Matlab docs. Creates a diagonal matrix with half its entries being negative.
            int n = evenOrder; 
            var A = Matrix.CreateZero(n, n);
            for (int i = 0; i < n / 2; ++i) A[i, i] = n / 2 - i;
            for (int i = 0; i < n / 2; ++i) A[n / 2 + i, n / 2 + i] = -1 - i;

            // x = [1, 1, ... 1]
            var x = Vector.CreateWithValue(n, 1.0);

            // b[i] = A[i, i] 
            var b = Vector.CreateZero(n);
            for (int i = 0; i < n; ++i) b[i] = A[i, i];

            // The usual Jacobi preconditioner worsens convergence. Modifications are needed.
            var positiveDiagonal = new double[n];
            for (int i = 0; i < n; ++i) positiveDiagonal[i] = Math.Abs(A[i, i]);
            var M = new JacobiPreconditioner(positiveDiagonal);
            return (A, b, x, M);
        }
    }
}
