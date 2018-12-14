using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Exceptions;

//TODO: I am pretty sure that the CSC calculations can use the corresponding CSR implementations by transposing the matrix. 
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

        //TODO: use transposed CSR method
        public void Dcscgemm(bool transposeA, int numRowsA, int numColsB, int numColsA, double[] valuesA, int[] colOffsetsA,
            int[] rowIndicesA, double[] b, double[] c)
            => Dcsrgemm(!transposeA, numColsA, numColsB, numRowsA, valuesA, colOffsetsA, rowIndicesA, b, c);

        public void Dcscgemv(bool transposeA, int numRowsA, int numColsA, double[] valuesA, int[] colOffsetsA, int[] rowIndicesA,
            double[] x, int offsetX, double[] y, int offsetY)
            => Dcsrgemv(!transposeA, numColsA, numRowsA, valuesA, colOffsetsA, rowIndicesA, x, offsetX, y, offsetY);

        public void Dcsrgemm(bool transposeA, int numRowsA, int numColsB, int numColsA, double[] valuesA, int[] rowOffsetsA,
            int[] colIndicesA, double[] b, double[] c)
        {
            int ldB = transposeA ? numRowsA : numColsA;
            int ldC = transposeA ? numColsA : numRowsA;

            //TODO: This implementation uses level 2 BLAS. I should implement it from scratch as a lvl 3 BLAS, but is it worth 
            //      it without parallelism?
            for (int j = 0; j < numColsB; ++j)
            {
                Dcsrgemv(transposeA, numRowsA, numColsA, valuesA, rowOffsetsA, colIndicesA, b, j * ldB, c, j * ldC);
            }
        }

        public void Dcsrgemv(bool transposeA, int numRowsA, int numColsA, double[] valuesA, int[] rowOffsetsA, int[] colIndicesA,
            double[] x, int offsetX, double[] y, int offsetY)
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

        /// <summary>
        /// Matrix vector multiplication y = A * x, with A being a symmetric matrix in Skyline format, where only the upper 
        /// triangle is stored.
        /// </summary>
        public void Dskymv(int order, double[] valuesA, int[] diagOffsetsA, double[] x, double[] y)
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
                double linearCombinationCoeff = x[j];
                // Dot product of the (L+D) part of the row * vector
                double dotLower = valuesA[diagOffset] * linearCombinationCoeff; // Contribution of diagonal entry: A[j,j] * x[j]
                for (int i = j - 1; i >= columnTop; --i) // Process the rest of the non zero entries of the column
                {
                    double aij = valuesA[diagOffset + j - i]; // Thus the matrix is only indexed once

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

        public void Dskysm(int order, int numRhs, double[] valuesA, int[] diagOffsetsA, double[] b, double[] x)
        {
            //TODO: This implementation uses level 2 BLAS. I should implement it from scratch as a lvl 3 BLAS, but is it worth 
            //      it without parallelism?
            for (int j = 0; j < numRhs; ++j)
            {
                int offset = j * order;
                Dskysv(order, valuesA, diagOffsetsA, b, offset, x, offset);
            }
        }

        /// <summary>
        /// Linear system solution x = inv(A) * b, with A being with a symmetric matrix in Skyline format, where only the upper 
        /// triangle is stored.
        /// </summary>
        public void Dskysv(int order, double[] valuesA, int[] diagOffsetsA, double[] b, double[] x)
        {
            // Copied from Stavroulakis code.

            // Copy the y vector
            Array.Copy(b, x, order);

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
                        C += valuesA[KK] * x[k];
                    }
                    x[n] -= C;
                }
            }

            // Back substitution
            for (n = 0; n < order; n++) x[n] /= valuesA[diagOffsetsA[n]];

            n = order - 1;
            for (int l = 1; l < order; l++)
            {
                int KL = diagOffsetsA[n] + 1;
                int KU = diagOffsetsA[n + 1] - 1;
                if (KU >= KL)
                {
                    int k = n;
                    double xn = x[n];
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        x[k] -= valuesA[KK] * xn;
                    }
                }
                n--;
            }
        }

        //TODO: this is the same method as above, but without the offsets, to avoid the extra computations. Perhaps I can 
        //      write this method, such that all index variables are initialized with respect to the offset, instead of
        //      always adding it. Also other micro optimizations should be doable.
        /// <summary>
        /// Linear system solution x = inv(A) * b, with A being with a symmetric matrix in Skyline format, where only the upper 
        /// triangle is stored.
        /// </summary>
        public void Dskysv(int order, double[] valuesA, int[] diagOffsetsA, double[] b, int offsetB, double[] x, int offsetX)
        {
            // Copied from Stavroulakis code.

            // Copy the y vector
            Array.Copy(b, offsetB, x, offsetX, order);

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
                        C += valuesA[KK] * x[k];
                    }
                    x[offsetX + n] -= C;
                }
            }

            // Back substitution
            for (n = 0; n < order; n++) x[offsetX + n] /= valuesA[diagOffsetsA[n]];

            n = order - 1;
            for (int l = 1; l < order; l++)
            {
                int KL = diagOffsetsA[n] + 1;
                int KU = diagOffsetsA[n + 1] - 1;
                if (KU >= KL) //TODO: can't I avoid this check, by accessing the entries in a smarter way?
                {
                    int k = offsetX + n;
                    double xn = x[k];
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        x[k] -= valuesA[KK] * xn;
                    }
                }
                n--;
            }
        }
    }
}
