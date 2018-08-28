using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Factorizations
{
    /// <summary>
    /// The Bunch-Kaufman factorization of a symmetric matrix A is" A = transpose(U) * D * U. It is recommended for symmetric 
    /// indefinite matrices, since Cholesky factorization cannot be used and Bunch-Kaufman is more efficient than LU 
    /// factorization. Uses Intel MKL.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class BunchKaufmanFactorization : ITriangulation
    {
        public double CalcDeterminant()
        {
            throw new NotImplementedException();
        }

        public Vector SolveLinearSystem(Vector rhs)
        {
            throw new NotImplementedException();
        }
    }
}
