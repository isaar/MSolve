using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    /// <summary>
    /// Uses managed C# code (usually unoptimized) to performs BLAS operations.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ManagedSparseBlasProvider : ISparseBlasProvider
    {
        public void CscTimesVector(int numCols, double[] values, int[] colOffsets, int[] rowIndices, double[] x, double[] y)
        {
            // A * x = linear combination of columns of A, with the entries of x as coefficients
            for (int j = 0; j < numCols; ++j)
            {
                double scalar = x[j];
                int colStart = colOffsets[j]; //inclusive
                int colEnd = colOffsets[j + 1]; //exclusive
                for (int k = colStart; k < colEnd; ++k)
                {
                    y[rowIndices[k]] += scalar * values[k];
                }
            }
        }

        public void CscTransposeTimesVector(int numCols, double[] values, int[] colOffsets, int[] rowIndices, double[] x, 
            double[] y)
        {
            // A^T * x = sum{row of A^T * x} = sum{col of A * x}
            for (int j = 0; j < numCols; ++j)
            {
                double dot = 0.0;
                int colStart = colOffsets[j]; //inclusive
                int colEnd = colOffsets[j + 1]; //exclusive
                for (int k = colStart; k < colEnd; ++k)
                {
                    dot += values[k] * x[rowIndices[k]];
                }
                y[j] = dot;
            }
        }

        public void CsrTimesVector(int numRows, double[] values, int[] rowOffsets, int[] colIndices, double[] x, double[] y)
        {
            for (int i = 0; i < numRows; ++i)
            {
                double dot = 0.0;
                int rowStart = rowOffsets[i]; //inclusive
                int rowEnd = rowOffsets[i + 1]; //exclusive
                for (int k = rowStart; k < rowEnd; ++k)
                {
                    dot += values[k] * x[colIndices[k]];
                }
                y[i] = dot;
            }
        }

        public void CsrTransposeTimesVector(int numRows, double[] values, int[] rowOffsets, int[] colIndices, double[] x, 
            double[] y)
        {
            // A^T * x = linear combination of columns of A^T = rows of A, with the entries of x as coefficients
            for (int i = 0; i < numRows; ++i)
            {
                double scalar = x[i];
                int rowStart = rowOffsets[i]; //inclusive
                int rowEnd = rowOffsets[i + 1]; //exclusive
                for (int k = rowStart; k < rowEnd; ++k)
                {
                    y[colIndices[k]] += scalar * values[k];
                }
            }
        }

        public void SkylineTimesVector(int order, double[] values, int[] diagOffsets, double[] x, double[] y)
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
                double linearCombinationCoeff = x[j];
                // Dot product of the (L+D) part of the row * vector
                double dotLower = values[diagOffset] * linearCombinationCoeff; // Contribution of diagonal entry: A[j,j] * x[j]
                for (int i = j - 1; i >= columnTop; --i) // Process the rest of the non zero entries of the column
                {
                    double aij = values[diagOffset + j - i]; // Thus the matrix is only indexed once

                    // Contribution of the L part of the row, which is identical to the stored column j.
                    // Thus A[j,i]=A[i,j] and sum(A[i,j]*x[j]) = sum(A[i,j]*x[i])
                    dotLower += aij * x[i];

                    // Contribution of the U part of the column: result += coefficient * column j of U. This will update all rows
                    // [columnTop, j) of the result vector need to be updated to account for the current column j. 
                    y[i] += aij * linearCombinationCoeff;
                }
                // Column j alters rows [0,j) of the result vector, thus this should be the 1st time result[j] is written.
                Debug.Assert(y[j] == 0);
                y[j] = dotLower; // contribution of the (L+D) part of the row. 
            }
        }
    }
}
