using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: many of these may be able to be simplified to avoid much code duplication. An approach would be to have dedicated 1D 
//      vectors that are rows or columns of matrices. Then the CSR loops will operate only on them. The rest of the methods, 
//      will just need to provide the correct rows/columns.
//TODO: These should not be in the Commons namespace.
namespace ISAAR.MSolve.LinearAlgebra.Commons
{
    /// <summary>
    /// Implementations of multiplication operations with a matrix stored in CSR format.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class CsrMultiplications
    {
        internal static void CsrTimesMatrix(int numCsrRows, double[] csrValues, int[] csrRowOffsets, int[] csrColIndices,
            IMatrixView other, Matrix result)
        {
            for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time
            {
                for (int i = 0; i < numCsrRows; ++i)
                {
                    double dot = 0.0;
                    int csrRowStart = csrRowOffsets[i]; //inclusive
                    int csrRowEnd = csrRowOffsets[i + 1]; //exclusive
                    for (int k = csrRowStart; k < csrRowEnd; ++k) // csr.row[i] * other.col[c]
                    {
                        dot += csrValues[k] * other[csrColIndices[k], c];
                    }
                    result[i, c] = dot;
                }
            }
        }

        internal static void CsrTimesMatrixTrans(int numCsrRows, double[] csrValues, int[] csrRowOffsets, int[] csrColIndices,
            IMatrixView other, Matrix result)
        {
            for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time
            {
                for (int i = 0; i < numCsrRows; ++i)
                {
                    double dot = 0.0;
                    int csrRowStart = csrRowOffsets[i]; //inclusive
                    int csrRowEnd = csrRowOffsets[i + 1]; //exclusive
                    for (int k = csrRowStart; k < csrRowEnd; ++k) // csr.row[i] * other.row[c]
                    {
                        dot += csrValues[k] * other[c, csrColIndices[k]];
                    }
                    result[i, c] = dot;
                }
            }
        }

        internal static void CsrTimesVector(int numCsrRows, double[] csrValues, int[] csrRowOffsets, int[] csrColIndices,
            IVectorView lhs, double[] rhs)
        {
            for (int i = 0; i < numCsrRows; ++i)
            {
                double dot = 0.0;
                int rowStart = csrRowOffsets[i]; //inclusive
                int rowEnd = csrRowOffsets[i + 1]; //exclusive
                for (int k = rowStart; k < rowEnd; ++k) dot += csrValues[k] * lhs[csrColIndices[k]];
                rhs[i] = dot;
            }
        }

        internal static void CsrTimesVector(int numCsrRows, double[] csrValues, int[] csrRowOffsets, int[] csrColIndices,
            IVectorView lhs, IVector rhs)
        {
            for (int i = 0; i < numCsrRows; ++i)
            {
                double dot = 0.0;
                int rowStart = csrRowOffsets[i]; //inclusive
                int rowEnd = csrRowOffsets[i + 1]; //exclusive
                for (int k = rowStart; k < rowEnd; ++k) dot += csrValues[k] * lhs[csrColIndices[k]];
                rhs.Set(i, dot);
            }
        }

        internal static void CsrTransTimesMatrix(int numCsrRows, double[] csrValues, int[] csrRowOffsets, int[] csrColIndices,
            IMatrixView other, Matrix result)
        {
            for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time
            {
                // A^T * x = linear combination of columns of A^T = rows of A, with the entries of x as coefficients, 
                // where x is column c of the other matrix
                for (int i = 0; i < numCsrRows; ++i)
                {
                    double scalar = other[i, c];
                    int csrRowStart = csrRowOffsets[i]; //inclusive
                    int csrRowEnd = csrRowOffsets[i + 1]; //exclusive
                    for (int k = csrRowStart; k < csrRowEnd; ++k) // sum(other[i,c] * transpose(csr.col[c])) = sum(other[i,c] * csr.row[c])
                    {
                        result[csrColIndices[k], c] += scalar * csrValues[k];
                    }
                }
            }
        }

        internal static void CsrTransTimesMatrixTrans(int numCsrRows, double[] csrValues, int[] csrRowOffsets,
            int[] csrColIndices, IMatrixView other, Matrix result)
        {
            for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time
            {
                // A^T * x = linear combination of columns of A^T = rows of A, with the entries of x as coefficients, 
                // where x is column c of the other matrix
                for (int i = 0; i < numCsrRows; ++i)
                {
                    double scalar = other[c, i];
                    int csrRowStart = csrRowOffsets[i]; //inclusive
                    int csrRowEnd = csrRowOffsets[i + 1]; //exclusive
                    for (int k = csrRowStart; k < csrRowEnd; ++k) // sum(other[c,i] * transpose(csr.col[c])) = sum(other[c,i] * csr.row[c])
                    {
                        result[csrColIndices[k], c] += scalar * csrValues[k];
                    }
                }
            }
        }

        internal static void CsrTransTimesVector(int numCsrRows, double[] csrValues, int[] csrRowOffsets, int[] csrColIndices,
            IVectorView lhs, double[] rhs)
        {
            // A^T * x = linear combination of columns of A^T = rows of A, with the entries of x as coefficients
            for (int i = 0; i < numCsrRows; ++i)
            {
                double scalar = lhs[i];
                int rowStart = csrRowOffsets[i]; //inclusive
                int rowEnd = csrRowOffsets[i + 1]; //exclusive
                for (int k = rowStart; k < rowEnd; ++k)
                {
                    rhs[csrColIndices[k]] += scalar * csrValues[k];
                }
            }
        }

        internal static void CsrTransTimesVector(int numCsrRows, double[] csrValues, int[] csrRowOffsets, int[] csrColIndices,
            IVectorView lhs, IVector rhs)
        {
            var temp = new double[rhs.Length];
            CsrTransTimesVector(numCsrRows, csrValues, csrRowOffsets, csrColIndices, lhs, temp);
            rhs.CopyFrom(Vector.CreateFromArray(temp));

            // The following requires a lot of indexing into the rhs vector.
            //// A^T * x = linear combination of columns of A^T = rows of A, with the entries of x as coefficients
            //for (int i = 0; i < numCsrRows; ++i)
            //{
            //    double scalar = lhs[i];
            //    int rowStart = csrRowOffsets[i]; //inclusive
            //    int rowEnd = csrRowOffsets[i + 1]; //exclusive
            //    for (int k = rowStart; k < rowEnd; ++k)
            //    {
            //        rhs.Set(csrColIndices[k], rhs[csrColIndices[k]] + scalar * csrValues[k]);
            //    }
            //}
        }

        internal static void MatrixTimesCsr(int numCsrRows, double[] csrValues, int[] csrRowOffsets, int[] csrColIndices,
            IMatrixView other, Matrix result)
        {
            for (int r = 0; r < result.NumRows; ++r) // Compute one output row at a time.
            {
                // x * A = linear combination of rows of A with the entries of x as coefficients,
                // where x is row r of the other matrix.
                for (int i = 0; i < numCsrRows; ++i)
                {
                    double scalar = other[r, i];
                    int csrRowStart = csrRowOffsets[i]; //inclusive
                    int csrRowEnd = csrRowOffsets[i + 1]; //exclusive
                    for (int k = csrRowStart; k < csrRowEnd; ++k) // sum(other[r,i] * csr.row[i]))
                    {
                        result[r, csrColIndices[k]] += scalar * csrValues[k];
                    }
                }
            }
        }

        internal static void MatrixTimesCsrTrans(int numCsrRows, double[] csrValues, int[] csrRowOffsets, int[] csrColIndices,
            IMatrixView other, Matrix result)
        {
            for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time.
            {
                int csrRowStart = csrRowOffsets[c]; //inclusive
                int csrRowEnd = csrRowOffsets[c + 1]; //exclusive
                for (int i = 0; i < other.NumRows; ++i)
                {
                    double dot = 0.0;
                    for (int k = csrRowStart; k < csrRowEnd; ++k) // other.row[i] * transpose(csr).col[j] = other.row[i] * csr.row[j]
                    {
                        dot += csrValues[k] * other[i, csrColIndices[k]];
                    }
                    result[i, c] = dot;
                }
            }
        }

        internal static void MatrixTransTimesCsr(int numCsrRows, double[] csrValues, int[] csrRowOffsets, int[] csrColIndices,
            IMatrixView other, Matrix result)
        {
            for (int r = 0; r < result.NumRows; ++r) // Compute one output row at a time.
            {
                // x * A = linear combination of rows of A with the entries of x as coefficients,
                // where x is row r of transpose(other matrix).
                for (int i = 0; i < numCsrRows; ++i)
                {
                    double scalar = other[i, r];
                    int csrRowStart = csrRowOffsets[i]; //inclusive
                    int csrRowEnd = csrRowOffsets[i + 1]; //exclusive
                    for (int k = csrRowStart; k < csrRowEnd; ++k) // sum(other[i,r] * csr.row[i]))
                    {
                        result[r, csrColIndices[k]] += scalar * csrValues[k];
                    }
                }
            }
        }

        internal static void MatrixTransTimesCsrTrans(int numCsrRows, double[] csrValues, int[] csrRowOffsets,
            int[] csrColIndices, IMatrixView other, Matrix result)
        {
            for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time
            {
                int rowStart = csrRowOffsets[c]; //inclusive
                int rowEnd = csrRowOffsets[c + 1]; //exclusive
                for (int i = 0; i < other.NumColumns; ++i)
                {
                    double dot = 0.0;
                    for (int k = rowStart; k < rowEnd; ++k) // other.col[i] * transpose(csr).col[j] = other.col[i] * csr.row[j]
                    {
                        dot += csrValues[k] * other[csrColIndices[k], i];
                    }
                    result[i, c] = dot;
                }
            }
        }
    }
}
