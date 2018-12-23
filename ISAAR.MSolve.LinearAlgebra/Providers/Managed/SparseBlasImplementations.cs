using System;
using System.Diagnostics;

namespace ISAAR.MSolve.LinearAlgebra.Providers.Managed
{
    /// <summary>
    /// Custom and unoptimized managed implementations of SparseBLAS operations, for which I have not found 3rd party 
    /// implementations yet. If performance is of concern, then the user should choose the native optimized providers (e.g. MKL).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class SparseBlasImplementations
    {
        internal static void AlphaTimesSparsePlusDenseVector(int nnz, double alpha, double[] x, int[] indicesX, int offsetX, 
            double[] y, int offsetY)
        {
            for (int i = 0; i < nnz; ++i)
            {
                int idxX = offsetX + i;
                y[offsetY + indicesX[idxX]] += alpha * x[idxX];
            }
        }

        internal static void CsrTimesVector(bool transposeA, int numRowsA, int numColsA, double[] valuesA, int[] rowOffsetsA, 
            int[] colIndicesA, double[] x, int offsetX, double[] y, int offsetY)
        {
            if (transposeA)
            {
                // A^T * x = linear combination of columns of A^T = rows of A, with the entries of x as coefficients
                for (int i = 0; i < numRowsA; ++i)
                {
                    double scalar = x[offsetX + i];
                    int rowStart = rowOffsetsA[i]; //inclusive
                    int rowEnd = rowOffsetsA[i + 1]; //exclusive
                    for (int k = rowStart; k < rowEnd; ++k)
                    {
                        y[offsetY + colIndicesA[k]] += scalar * valuesA[k];
                    }
                }
            }
            else
            {
                for (int i = 0; i < numRowsA; ++i)
                {
                    double dot = 0.0;
                    int rowStart = rowOffsetsA[i]; //inclusive
                    int rowEnd = rowOffsetsA[i + 1]; //exclusive
                    for (int k = rowStart; k < rowEnd; ++k) dot += valuesA[k] * x[offsetX + colIndicesA[k]];
                    y[offsetY + i] = dot;
                }
            }
        }

        internal static void SkylineTimesVector(int order, double[] valuesA, int[] diagOffsetsA, double[] vectorX, 
            double[] vectorY)
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
                int diagOffset = diagOffsetsA[j];
                int columnTop = j - diagOffsetsA[j + 1] + diagOffset + 1;
                double linearCombinationCoeff = vectorX[j];
                // Dot product of the (L+D) part of the row * vector
                double dotLower = valuesA[diagOffset] * linearCombinationCoeff; // Contribution of diagonal entry: A[j,j] * x[j]
                for (int i = j - 1; i >= columnTop; --i) // Process the rest of the non zero entries of the column
                {
                    double aij = valuesA[diagOffset + j - i]; // Thus the matrix is only indexed once

                    // Contribution of the L part of the row, which is identical to the stored column j.
                    // Thus A[j,i]=A[i,j] and sum(A[i,j]*x[j]) = sum(A[i,j]*x[i])
                    dotLower += aij * vectorX[i];

                    // Contribution of the U part of the column: result += coefficient * column j of U. This will update all rows
                    // [columnTop, j) of the result vector need to be updated to account for the current column j. 
                    vectorY[i] += aij * linearCombinationCoeff;
                }
                // Column j alters rows [0,j) of the result vector, thus this should be the 1st time result[j] is written.
                Debug.Assert(vectorY[j] == 0);
                vectorY[j] = dotLower; // contribution of the (L+D) part of the row. 
            }
        }

        internal static void SkylineSystemSolution(int order, double[] valuesA, int[] diagOffsetsA, double[] vectorB, 
            double[] vectorX)
        {
            // Copied from Stavroulakis code.

            // Copy the y vector
            Array.Copy(vectorB, vectorX, order);

            // RHS vector reduction
            int n;
            for (n = 0; n < order; n++)
            {
                int KL = diagOffsetsA[n] + 1;
                int KU = diagOffsetsA[n + 1] - 1;
                if (KU >= KL)
                {
                    int k = n;
                    double C = 0;
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        C += valuesA[KK] * vectorX[k];
                    }
                    vectorX[n] -= C;
                }
            }

            // Back substitution
            for (n = 0; n < order; n++) vectorX[n] /= valuesA[diagOffsetsA[n]];

            n = order - 1;
            for (int l = 1; l < order; l++)
            {
                int KL = diagOffsetsA[n] + 1;
                int KU = diagOffsetsA[n + 1] - 1;
                if (KU >= KL)
                {
                    int k = n;
                    double xn = vectorX[n];
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        vectorX[k] -= valuesA[KK] * xn;
                    }
                }
                n--;
            }
        }

        //TODO: this is the same method as above, but without the offsets, to avoid the extra computations. Perhaps I can 
        //      write this method, such that all index variables are initialized with respect to the offset, instead of
        //      always adding it. Also other micro optimizations should be doable.
        internal static void SkylineSystemSolutionWithOffsets(int order, double[] valuesA, int[] diagOffsetsA, 
            double[] vectorB, int offsetB, double[] vectorX, int offsetX)
        {
            // Copied from Stavroulakis code.

            // Copy the y vector
            Array.Copy(vectorB, offsetB, vectorX, offsetX, order);

            // RHS vector reduction
            int n;
            for (n = 0; n < order; n++)
            {
                int KL = diagOffsetsA[n] + 1;
                int KU = diagOffsetsA[n + 1] - 1;
                if (KU >= KL) //TODO: can't I avoid this check, by accessing the entries in a smarter way?
                {
                    int k = offsetX + n;
                    double C = 0;
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        C += valuesA[KK] * vectorX[k];
                    }
                    vectorX[offsetX + n] -= C;
                }
            }

            // Back substitution
            for (n = 0; n < order; n++) vectorX[offsetX + n] /= valuesA[diagOffsetsA[n]];

            n = order - 1;
            for (int l = 1; l < order; l++)
            {
                int KL = diagOffsetsA[n] + 1;
                int KU = diagOffsetsA[n + 1] - 1;
                if (KU >= KL) //TODO: can't I avoid this check, by accessing the entries in a smarter way?
                {
                    int k = offsetX + n;
                    double xn = vectorX[k];
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        vectorX[k] -= valuesA[KK] * xn;
                    }
                }
                n--;
            }
        }

        internal static double SparseDotDenseVector(int nnz, double[] x, int[] indicesX, int offsetX, double[] y, int offsetY)
        {
            double dot = 0.0;
            for (int i = 0; i < nnz; ++i)
            {
                int idxX = offsetX + i;
                dot += x[idxX] * y[offsetY + indicesX[idxX]];
            }
            return dot;
        }
    }
}
