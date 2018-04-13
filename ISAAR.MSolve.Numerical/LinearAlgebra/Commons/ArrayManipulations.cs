using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Commons
{
    internal static class ArrayManipulations
    {
        /// <summary>
        /// Add zero columns to the right of the matrix or remove existing ones.
        /// </summary>
        /// <param name="numRows">The number of rows in <paramref name="matrix"/>.</param>
        /// <param name="numColsOld">The number of columns in <paramref name="matrix"/>.</param>
        /// <param name="numColsNew">The number of columns after the operation.</param>
        /// <param name="matrix">The matrix in full column major layout.</param>
        /// <returns></returns>
        internal static double[] ResizeColsColMajor(int numRows, int numColsOld, int numColsNew, double[] matrix)
        {
            // We need to copy the common columns. Since the matrix is column major, they are contiguous
            double[] result = new double[numRows * numColsNew];
            int numColsMin = (numColsNew > numColsOld) ? numColsOld : numColsNew;
            Array.Copy(matrix, result, matrix.Length);
            return result;
        }

        /// <summary>
        /// Add zero rows to the bottom of the matrix.
        /// </summary>
        /// <param name="numRowsOld">The number of rows in <paramref name="matrix"/>.</param>
        /// <param name="numCols">The number of columns in <paramref name="matrix"/></param>
        /// <param name="numRowsNew">The number of columns after the operation. It must be 
        ///     <paramref name="numRowsNew"/> &gt; <paramref name="numRowsOld"/>.</param>
        /// <param name="matrix">The matrix in full column major layout.</param>
        /// <returns></returns>
        internal static double[] IncreaseRowsColMajor(int numRowsOld, int numCols, int numRowsNew, double[] matrix)
        {
            // We need to copy the common rows. Since the matrix is column major, they will not be contiguous.
            // Instead copy one column of length = numRowsOld at a time. Its entries are contiguous.
            double[] result = new double[numRowsNew * numCols];
            int contiguousLength = numRowsOld;
            for (int j = 0; j < numCols; ++j) 
            {
                int sourceStart = j * numRowsOld;
                int destinationStart = j * numRowsNew;
                Array.Copy(matrix, sourceStart, result, destinationStart, contiguousLength);
            }
            return result;
        }
    }
}
