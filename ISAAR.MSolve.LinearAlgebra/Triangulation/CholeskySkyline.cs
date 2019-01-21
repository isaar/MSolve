using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Triangulation
{
    public class CholeskySkyline
    {
        // TODO: fix this. GetFactorU() assumes an A = L^T*D*L factorization, but this class uses A=L^T*L (I think).
        /// <summary>
        /// Explicitly creates the upper triangular matrix U that resulted from the Cholesky factorization: A = transpose(U) * U,
        /// where A and U are n-by-n. 
        /// This method is safe to use as the factorization data are copied (if necessary). However, it is inefficient if the 
        /// generated matrix is only used once.
        /// </summary>
        //public TriangularUpper GetFactorU()
        //{
        //    // The factorization A = transpose(u) * D * u, u = unit upper triangular is stored. Thus U = sqrt(D) * u.
        //    // Since D is diagonal, we need to scale each column j of u by sqrt(D[j,j]).
        //    var upper = TriangularUpper.CreateZero(NumColumns); 
        //    for (int j = 0; j < NumColumns; ++j)
        //    {
        //        int diagOffset = diagOffsets[j];
        //        double sqrtD = Math.Sqrt(values[diagOffset]); // The diagonal of D is stored instead of the unit diagonal entries of u
        //        int colTop = j - diagOffsets[j + 1] + diagOffset + 1;
        //        upper[j, j] = sqrtD; // U[j,j] = u[j,j] * sqrt(D[j,j]) = 1 * sqrt(D[j,j])
        //        for (int i = j - 1; i >= colTop; --i)
        //        {
        //            int offset = diagOffset + j - i;
        //            upper[i, j] = sqrtD * values[offset]; // U[:,j] = sqrt(D[j,j]) * u[:,j]
        //        }
        //    }
        //    int[] diagOffsetsCopy = new int[diagOffsets.Length];
        //    Array.Copy(diagOffsets, diagOffsetsCopy, diagOffsets.Length);
        //    return upper;
        //}
    }
}
