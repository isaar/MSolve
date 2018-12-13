using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Exceptions;

//TODO: perhaps custom implementations should be in a dedicated namespace. Or they shoudl be in the providers. Keep one design 
//      though.
//TODO: At some point I should change my Skyline format to match the one used in MKL. Then the SparseBLAS operations can be 
//      interchangeable.
//TODO: This should be a singleton
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    /// <summary>
    /// Uses managed C# code (usually unoptimized) to perform BLAS operations.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ManagedSparseBlasProvider : ISparseBlasProvider
    {
        public void Daxpyi(int nnz, double[] alpha, double[] x, int offsetX, int[] indicesX, double[] y, int offsetY)
        {
            throw new NotImplementedException();
        }

        public void Dcscgemv(bool transpose, int numRows, int numCols, double[] values, int[] colOffsets, int[] rowIndices,
            double[] lhs, double[] rhs)
        {
            if (transpose)
            {
                // A^T * x = sum{row of A^T * x} = sum{col of A * x}
                for (int j = 0; j < numCols; ++j)
                {
                    double dot = 0.0;
                    int colStart = colOffsets[j]; //inclusive
                    int colEnd = colOffsets[j + 1]; //exclusive
                    for (int k = colStart; k < colEnd; ++k)
                    {
                        dot += values[k] * lhs[rowIndices[k]];
                    }
                    rhs[j] = dot;
                }
            }
            else
            {
                // A * x = linear combination of columns of A, with the entries of x as coefficients
                for (int j = 0; j < numCols; ++j)
                {
                    double scalar = lhs[j];
                    int colStart = colOffsets[j]; //inclusive
                    int colEnd = colOffsets[j + 1]; //exclusive
                    for (int k = colStart; k < colEnd; ++k)
                    {
                        rhs[rowIndices[k]] += scalar * values[k];
                    }
                }
            }
        }

        public void Dcsrgemv(bool transpose, int numRows, int numColumns, double[] values, int[] rowOffsets, int[] colIndices,
            double[] lhs, double[] rhs)
        {
            if (transpose)
            {
                // A^T * x = linear combination of columns of A^T = rows of A, with the entries of x as coefficients
                for (int i = 0; i < numRows; ++i)
                {
                    double scalar = lhs[i];
                    int rowStart = rowOffsets[i]; //inclusive
                    int rowEnd = rowOffsets[i + 1]; //exclusive
                    for (int k = rowStart; k < rowEnd; ++k)
                    {
                        rhs[colIndices[k]] += scalar * values[k];
                    }
                }
            }
            else
            {
                for (int i = 0; i < numRows; ++i)
                {
                    double dot = 0.0;
                    int rowStart = rowOffsets[i]; //inclusive
                    int rowEnd = rowOffsets[i + 1]; //exclusive
                    for (int k = rowStart; k < rowEnd; ++k)
                    {
                        dot += values[k] * lhs[colIndices[k]];
                    }
                    rhs[i] = dot;
                }
            }
        }

        /// <summary>
        /// Matrix vector multiplication, with a symmetric matrix in Skyline format, where only the upper triangle is stored.
        /// </summary>
        public void Dskymv(int order, double[] values, int[] diagOffsets, double[] lhs, double[] rhs)
        {
            // A*x = (L+D)*x + U*x
            // (L+D)*x is easy, since the non zero entries of row i left of the diagonal are stored contiguously in column i and
            // we can easily take its dot product with the vector.
            // U*x is trickier, since we cannot access contiguously the non zero entries of row i. Instead think of it as
            // U*x = linear combination of columns of U (accessed contiguously) with the entries of vector as coefficients. Then 
            // we can deal with them while we process the next columns (i, n-1]. This way the matrix is only indexed once, but 
            // not the result vector entry result[i].
            for (int j = 0; j < order; ++j)
            {
                int diagOffset = diagOffsets[j];
                int columnTop = j - diagOffsets[j + 1] + diagOffset + 1;
                double linearCombinationCoeff = lhs[j];
                // Dot product of the (L+D) part of the row * vector
                double dotLower = values[diagOffset] * linearCombinationCoeff; // Contribution of diagonal entry: A[j,j] * x[j]
                for (int i = j - 1; i >= columnTop; --i) // Process the rest of the non zero entries of the column
                {
                    double aij = values[diagOffset + j - i]; // Thus the matrix is only indexed once

                    // Contribution of the L part of the row, which is identical to the stored column j.
                    // Thus A[j,i]=A[i,j] and sum(A[i,j]*x[j]) = sum(A[i,j]*x[i])
                    dotLower += aij * lhs[i];

                    // Contribution of the U part of the column: result += coefficient * column j of U. This will update all rows
                    // [columnTop, j) of the result vector need to be updated to account for the current column j. 
                    rhs[i] += aij * linearCombinationCoeff;
                }
                // Column j alters rows [0,j) of the result vector, thus this should be the 1st time result[j] is written.
                Debug.Assert(rhs[j] == 0);
                rhs[j] = dotLower; // contribution of the (L+D) part of the row. 
            }
        }

        /// <summary>
        /// Linear system solution, with a symmetric matrix in Skyline format, where only the upper triangle is stored.
        /// </summary>
        public void Dskysv(int order, double[] values, int[] diagOffsets, double[] rhs, double[] lhs)
        {
            // Copied from Stavroulakis code.

            // Copy the rhs vector
            Array.Copy(rhs, lhs, order);

            // RHS vector reduction
            int n;
            for (n = 0; n < order; n++)
            {
                int KL = diagOffsets[n] + 1;
                int KU = diagOffsets[n + 1] - 1;
                if (KU >= KL)
                {
                    int k = n;
                    double C = 0;
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        C += values[KK] * lhs[k];
                    }
                    lhs[n] -= C;
                }
            }

            // Back substitution
            for (n = 0; n < order; n++) lhs[n] /= values[diagOffsets[n]];

            n = order - 1;
            for (int l = 1; l < order; l++)
            {
                int KL = diagOffsets[n] + 1;
                int KU = diagOffsets[n + 1] - 1;
                if (KU >= KL)
                {
                    int k = n;
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        lhs[k] -= values[KK] * lhs[n];
                    }
                }
                n--;
            }
        }
    }
}
