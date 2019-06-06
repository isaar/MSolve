using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Can't I just do an LDL for symmetric indefinite matrices? It works and it is as efficient as Cholesky.
namespace ISAAR.MSolve.LinearAlgebra.Triangulation
{
    /// <summary>
    /// The Bunch-Kaufman factorization of a symmetric matrix A is" A = transpose(U) * D * U. It is recommended for symmetric 
    /// indefinite matrices, since Cholesky factorization cannot be used and Bunch-Kaufman is more efficient than LU 
    /// factorization. Uses LAPACK.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class BunchKaufmanFactorization : ITriangulation
    {
        public int Order => throw new NotImplementedException();

        public double CalcDeterminant()
        {
            throw new NotImplementedException();
        }

        public void SolveLinearSystem(Vector rhs, Vector solution)
        {
            throw new NotImplementedException();
        }
    }
}
