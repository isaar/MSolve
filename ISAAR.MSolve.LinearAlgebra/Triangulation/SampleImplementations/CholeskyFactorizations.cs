using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.LinearAlgebra.Triangulation.SampleImplementations
{
    internal static class CholeskyFactorizations
    {
        /// <summary>
        /// Returns 0 if the factorization was successful. If the matrix is not positive definite, returns the column index 
        /// where a negative subroot appears first.
        /// </summary>
        internal static int FullUpper1(int order, Matrix matrix, double pivotTolerance) 
        {
            // I think this is called jik dot form
            //TODO: Efficiency can be improved by incrementing the offsets and using temp variables for some array reads.

            // Process column j
            for (int j = 0; j < order; ++j)
            {
                // Update each A[i,j] by traversing the column j from top until the diagonal.
                for (int i = 0; i < j; ++i)
                {
                    // Dot product of the parts of columns j, i (i < j) between: [row 0, row i) 
                    double dotColsIJ = 0.0;
                    for (int k = 0; k < i; ++k)
                    {
                        dotColsIJ += matrix[k, i] * matrix[k, j];
                    }
                    matrix[i, j] = (matrix[i, j] - dotColsIJ) / matrix[i, i];
                }

                // Update the diagonal term
                // Dot product with itself of the part of column j: between: [row 0, row j) 
                double dotColsJJ = 0.0;
                for (int k = 0; k < j; ++k)
                {
                    dotColsJJ += Math.Pow(matrix[k, j], 2);
                }

                // if A[j,j] = sqrt(A[j,j]-dotColsJJ), but if the subroot is <= 0, then the matrix is not positive definite
                double subroot = matrix[j, j] - dotColsJJ;
                if (subroot < pivotTolerance) return j;
                matrix[j, j] = Math.Sqrt(subroot);
            }
            return 0;
        }

        /// <summary>
        /// Returns 0 if the factorization was successful. If the matrix is not positive definite, returns the column index 
        /// where a negative subroot appears first.
        /// </summary>
        internal static int FullUpper2(int order, Matrix matrix, double pivotTolerance)
        {
            // This version can be found in 
            // https://algowiki-project.org/en/Cholesky_decomposition#Software_implementation_of_the_algorithm.
            // Since it traverses each row of U instead of its columns, it cannot be used with skyline format.
            //TODO: Efficiency can be improved by incrementing the offsets and using temp variables for some array reads.

            // Process row i
            for (int i = 0; i < order; ++i)
            {
                // Update the diagonal term
                double dotColsII = 0;
                for (int k = 0; k < i; ++k)
                {
                    dotColsII += Math.Pow(matrix[k, i], 2);
                }

                // if A[i,i] = sqrt(A[i,i]-dotII), but if the subroot is <= 0, then the matrix is not positive definite
                double subroot = matrix[i, i] - dotColsII;
                if (subroot < pivotTolerance) return i;
                matrix[i, i] = Math.Sqrt(subroot);

                // Update each A[i,j] by traversing the row i from the diagonal to the right
                for (int j = i+1; j < order; ++j)
                {
                    double dotColsIJ = 0.0;
                    for (int k = 0; k < i; ++k)
                    {
                        dotColsIJ += matrix[k, i] * matrix[k, j];
                    }
                    // A[i,j] = (A[i,j] - dotColsIJ) / A[i,i]
                    matrix[i, j] = (matrix[i, j] - dotColsIJ) / matrix[i, i];
                }
            }
            return 0;
        }

        internal static void SolveSystemFullUpper(int order, Matrix upper, double[] rhs, double[] sol)
        {
            // A * x = b => U^T * U * x = b => U^T * y = b and U * x = y
            Array.Copy(rhs, sol, rhs.Length);

            // Forward substitution: U^T * y = b
            for (int i = 0; i < order; ++i)
            {
                double dot = 0;
                for (int j = 0; j < i; ++j)
                {
                    // Since upper is column major, access to entries of column i is contiguous
                    dot += upper[j, i] * sol[j];
                }
                sol[i] = (sol[i] - dot) / upper[i, i];
            }

            // Back substitution
            for (int i = order - 1; i >= 0; --i)
            {
                double dot = 0;
                for (int j = i + 1; j < order; ++j)
                {
                    // Since upper is column major, access to entries of column i is NOT contiguous
                    dot += upper[i, j] * sol[j];
                }
                sol[i] = (sol[i] - dot) / upper[i, i];
            }
        }
    }
}
