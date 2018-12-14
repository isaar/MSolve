using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: most of these can be simplified.
namespace ISAAR.MSolve.LinearAlgebra.Commons
{
    /// <summary>
    /// Implementations of common BLAS like operations with a matrix stored in CSR format.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class CsrCscStrategies
    {


        internal static void CscTimesMatrix(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
            IMatrixView other, Matrix result)
        {
            for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time.
            {
                // A * x = linear combination of columns of A, with the entries of x as coefficients, 
                // where x is column c of the other matrix
                for (int j = 0; j < numCscCols; ++j)
                {
                    double scalar = other[j, c];
                    int cscColStart = cscColOffsets[j]; //inclusive
                    int cscColEnd = cscColOffsets[j + 1]; //exclusive
                    for (int k = cscColStart; k < cscColEnd; ++k) // sum(other[j,c] * csc.col[j]))
                    {
                        result[cscRowIndices[k], c] += scalar * cscValues[k];
                    }
                }
            }
        }

        internal static void CscTimesMatrixTrans(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
            IMatrixView other, Matrix result)
        {
            for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time.
            {
                // A * x = linear combination of columns of A, with the entries of x as coefficients, 
                // where x is column c of transpose(other matrix)
                for (int j = 0; j < numCscCols; ++j)
                {
                    double scalar = other[c, j];
                    int cscColStart = cscColOffsets[j]; //inclusive
                    int cscColEnd = cscColOffsets[j + 1]; //exclusive
                    for (int k = cscColStart; k < cscColEnd; ++k) // sum(other[c,j] * csc.col[j]))
                    {
                        result[cscRowIndices[k], c] += scalar * cscValues[k];
                    }
                }
            }
        }

        internal static void CscTransTimesMatrix(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
            IMatrixView other, Matrix result)
        {
            for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time.
            {
                for (int j = 0; j < numCscCols; ++j)
                {
                    int cscColStart = cscColOffsets[j]; //inclusive
                    int cscColEnd = cscColOffsets[j + 1]; //exclusive
                    double dot = 0.0;
                    for (int k = cscColStart; k < cscColEnd; ++k) // traspose(csc).row[i] * other.col[c]
                    {
                        dot += cscValues[k] * other[cscRowIndices[k], c];
                    }
                    result[j, c] = dot;
                }
            }
        }

        internal static void CscTransTimesMatrixTrans(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
            IMatrixView other, Matrix result)
        {
            for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time.
            {
                for (int j = 0; j < numCscCols; ++j)
                {
                    int cscColStart = cscColOffsets[j]; //inclusive
                    int cscColEnd = cscColOffsets[j + 1]; //exclusive
                    double dot = 0.0;
                    for (int k = cscColStart; k < cscColEnd; ++k) // traspose(csc).row[i] * other.col[c]
                    {
                        dot += cscValues[k] * other[c, cscRowIndices[k]];
                    }
                    result[j, c] = dot;
                }
            }
        }

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
    }
}
