using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra
{
    /// <summary>
    /// Defines conversions from vectors, matrices, single- and multi-dimensional double arrays to integer arrays of the same 
    /// dimensions and vice-versa.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class IntegerArrayExtensions
    {
        /// <summary>
        /// Copies an 1D double[] array to an 1D int[] array. Only the integer part of each double entry is retained.
        /// </summary>
        public static int[] ToIntArray(this double[] array)
        {
            var result = new int[array.Length];
            for (int i = 0; i < array.Length; ++i) result[i] = (int)array[i];
            return result;
        }

        /// <summary>
        /// Copies a 2D double[,] array to a 2D int[,] array. Only the integer part of each double entry is retained.
        /// </summary>
        public static int[,] ToIntArray(this double[,] array)
        {
            var result = new int[array.GetLength(0), array.GetLength(1)];
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                for (int j = 0; j < array.GetLength(1); ++j) result[i, j] = (int)array[i, j];
            }
            return result;
        }

        /// <summary>
        /// Copies a 3D double[,,] array to a 3D int[,,] array. Only the integer part of each double entry is retained.
        /// </summary>
        public static int[,,] ToIntArray(this double[,,] array)
        {
            var result = new int[array.GetLength(0), array.GetLength(1), array.GetLength(2)];
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                for (int j = 0; j < array.GetLength(1); ++j)
                {
                    for (int k = 0; k < array.GetLength(2); ++k)  result[i, j, k] = (int)array[i, j, k];
                }
            }
            return result;
        }

        /// <summary>
        /// Copies a vector to an 1D dimensional int[] array. Only the integer part of each double entry is retained.
        /// </summary>
        public static int[] ToIntArray(this IVectorView vector) //TODO: Also for IIndexable1D
            => vector.CopyToArray().ToIntArray();

        /// <summary>
        /// Copies a matrix to a 2D int[,] array. Only the integer part of each double entry is retained.
        /// </summary>
        public static int[,] ToIntArray(this IMatrixView matrix) //TODO: Also for IIndexable1D
            => matrix.CopytoArray2D().ToIntArray();
    }
}
