using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Calculation of various matrix norms.
    /// If one of the following methods is also implemented in a concrete class, use that one, as it will be far more efficient.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class MatrixNormExtensions
    {
        /// <summary>
        /// Calculates the infinite induced norm of a matrix, namely the maximum absolute row sum. See 
        /// https://www.mathworks.com/help/matlab/ref/norm.html#bvhji30-4
        /// </summary>
        /// <param name="matrix"></param>
        public static double NormInf(this IMatrixView matrix) //TODO: matrix norms should be in a dedicated extensions class.
        {
            double max = double.MinValue;
            for (int i = 0; i < matrix.NumRows; ++i)
            {
                double rowSum = 0;
                for (int j = 0; j < matrix.NumColumns; ++j) rowSum += Math.Abs(matrix[i, j]);
                if (rowSum > max) max = rowSum;
            }
            return max;
        }
    }
}
