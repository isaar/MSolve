using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.LinearAlgebra.Triangulation.SampleImplementations
{
    /// <summary>
    /// Basic method implementations for a direct solver using Cholesky factorization. The matrix is assumed to be stored in 
    /// some column major format and we work with the upper factor. 
    /// To solve a linear system A * x = b:  
    /// 1) Factorize(): A = U^T * U 
    /// 2) ForwardSubstitution(): U^T * y = b 
    /// 3) BackSubstituttion(): U * x = y
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class CholeskyFactorizations
    {
        /// <summary>
        /// Returns 0 if the factorization was successful. If the matrix is not positive definite, returns the column index 
        /// where a negative subroot appears first.
        /// </summary>
        internal static int FactorizeFullUpper1(int n, Matrix A, double pivotTolerance) 
        {
            // I think this is called jik dot form
            //TODO: Efficiency can be improved by incrementing the offsets and using temp variables for some array reads.

            // Process column j
            for (int j = 0; j < n; ++j)
            {
                // Update each A[i,j] by traversing the column j from top until the diagonal.
                for (int i = 0; i < j; ++i)
                {
                    // Dot product of the parts of columns j, i (i < j) between: [row 0, row i) 
                    double dotColsIJ = 0.0;
                    for (int k = 0; k < i; ++k)
                    {
                        dotColsIJ += A[k, i] * A[k, j];
                    }
                    A[i, j] = (A[i, j] - dotColsIJ) / A[i, i];
                }

                // Update the diagonal term
                // Dot product with itself of the part of column j: between: [row 0, row j) 
                double dotColsJJ = 0.0;
                for (int k = 0; k < j; ++k)
                {
                    dotColsJJ += Math.Pow(A[k, j], 2);
                }

                // if A[j,j] = sqrt(A[j,j]-dotColsJJ), but if the subroot is <= 0, then the matrix is not positive definite
                double subroot = A[j, j] - dotColsJJ;
                if (subroot < pivotTolerance) return j;
                A[j, j] = Math.Sqrt(subroot);
            }
            return 0;
        }

        /// <summary>
        /// Returns 0 if the factorization was successful. If the matrix is not positive definite, returns the column index 
        /// where a negative subroot appears first.
        /// </summary>
        internal static int FactorizeFullUpper2(int n, Matrix A, double pivotTolerance)
        {
            // This version can be found in 
            // https://algowiki-project.org/en/Cholesky_decomposition#Software_implementation_of_the_algorithm.
            // Since it traverses each row of U instead of its columns, it cannot be used with skyline format.
            //TODO: Efficiency can be improved by incrementing the offsets and using temp variables for some array reads.

            // Process row i
            for (int i = 0; i < n; ++i)
            {
                // Update the diagonal term
                double dotColsII = 0;
                for (int k = 0; k < i; ++k)
                {
                    dotColsII += Math.Pow(A[k, i], 2);
                }

                // if A[i,i] = sqrt(A[i,i]-dotII), but if the subroot is <= 0, then the matrix is not positive definite
                double subroot = A[i, i] - dotColsII;
                if (subroot < pivotTolerance) return i;
                A[i, i] = Math.Sqrt(subroot);

                // Update each A[i,j] by traversing the row i from the diagonal to the right
                for (int j = i+1; j < n; ++j)
                {
                    double dotColsIJ = 0.0;
                    for (int k = 0; k < i; ++k)
                    {
                        dotColsIJ += A[k, i] * A[k, j];
                    }
                    // A[i,j] = (A[i,j] - dotColsIJ) / A[i,i]
                    A[i, j] = (A[i, j] - dotColsIJ) / A[i, i];
                }
            }
            return 0;
        }

        internal static void ForwardSubstitutionFullUpper(int n, Matrix U, double[] b, double[] x)
        {
            for (int i = 0; i < n; ++i)
            {
                double dot = 0;
                for (int j = 0; j < i; ++j)
                {
                    // Since upper is column major, access to entries of column i is contiguous
                    dot += U[j, i] * x[j];
                }
                x[i] = (b[i] - dot) / U[i, i];
            }
        }

        internal static void BackSubstitutionFullUpper1(int n, Matrix U, double[] x)
        {
            // This is version of the algorithm could be called vector / column / axpy. Dongarra provides various names for
            // matrix-vector multiplications and this is similar, so search there. It can be extended to sparse matrices, where 
            // the upper triangle is stored in column major sparse formats (e.g. skyline).

            // In this version, we will process the matrix column-by-column (dot version does row-by-row). Once we calculate the
            // diagonal entry we will scale and subtract the reduced column (above the diagonal) from the remaining unknowns in x.

            for (int j = n - 1; j >= 0; --j)
            {
                x[j] /= U[j, j];
                double diagValue = -x[j];
                for (int i = 0; i < j; ++i)
                {
                    x[i] += diagValue * U[i, j];
                }
            }
        }

        internal static void BackSubstitutionFullUpper2(int n, Matrix U, double[] x)
        {
            // This is the dot version of the algorithm. It is straight forward and probably the most efficient. However it 
            // cannot be extended to sparse matrices, where the upper triangle is stored in column major sparse formats 
            // (e.g. skyline).

            for (int i = n - 1; i >= 0; --i)
            {
                double dot = 0;
                for (int j = i + 1; j < n; ++j)
                {
                    // Since upper is column major, access to entries of column i is NOT contiguous
                    dot += U[i, j] * x[j];
                }
                x[i] = (x[i] - dot) / U[i, i];
            }
        }
    }
}
