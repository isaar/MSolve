using System;
using System.Collections.Generic;
using System.Text;

//TODO: this and matrix extensions whould not be in the same folder with the actual matrices and their interfaces
namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    public static class SparseMatrixExtensions
    {
        /// <summary>
        /// Only entries that satisfy A[i,j] != 0.0 will be included.
        /// </summary>
        /// <param name="matrix"></param>
        /// <returns></returns>
        public static IEnumerable<(int row, int col, double val)> EnumerateNonZeros(this double[,] matrix)
        {
            for (int i = 0; i < matrix.GetLength(0); ++i)
            {
                for (int j = 0; j < matrix.GetLength(1); ++j)
                {
                    if (matrix[i, j] != 0.0) yield return (i, j, matrix[i, j]);
                }
            }
        }

        /// <summary>
        /// Only entries that satisfy Math.Abs(A[i, j]) &gt; tolerance will be included.
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        public static IEnumerable<(int row, int col, double val)> EnumerateNonZeros(this double[,] matrix, double tolerance)
        {
            for (int i = 0; i < matrix.GetLength(0); ++i)
            {
                for (int j = 0; j < matrix.GetLength(1); ++j)
                {
                    if (Math.Abs(matrix[i, j]) > 0.0) yield return (i, j, matrix[i, j]);
                }
            }
        }

        /// <summary>
        /// Only entries that satisfy A[i,j] != 0.0 will be included.
        /// </summary>
        /// <param name="matrix"></param>
        /// <returns></returns>
        public static IEnumerable<(int row, int col, double val)> EnumerateNonZeros(this Matrix matrix)
        {
            for (int j = 0; j < matrix.NumColumns; ++j)
            {
                for (int i = 0; i < matrix.NumRows; ++i)
                {
                    if (matrix[i, j] != 0.0) yield return (i, j, matrix[i, j]);
                }
            }
        }

        /// <summary>
        /// Only entries that satisfy Math.Abs(A[i, j]) &gt; tolerance will be included.
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        public static IEnumerable<(int row, int col, double val)> EnumerateNonZeros(this Matrix matrix, double tolerance)
        {
            for (int j = 0; j < matrix.NumColumns; ++j)
            {
                for (int i = 0; i < matrix.NumRows; ++i)
                {
                    if (Math.Abs(matrix[i, j]) > 0.0) yield return (i, j, matrix[i, j]);
                }
            }
        }
    }
}
