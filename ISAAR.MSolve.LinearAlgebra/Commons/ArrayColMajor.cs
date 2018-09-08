using System;
using System.Collections.Generic;

//TODO: These should be delegated to C .dlls or to MKL if possible.
namespace ISAAR.MSolve.LinearAlgebra.Commons
{
    /// <summary>
    /// Low level array operations for matrices stored in full column major layout.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class ArrayColMajor
    {
        /// <summary>
        /// result = [<paramref name="matrixLeft"/>, <paramref name="matrixRight"/>] (Matlab notation). Dimensions of result:
        /// (<paramref name="numRowsLeft"/>=<paramref name="numRowsRight"/>)-by-
        /// (<paramref name="numColsLeft"/>+<paramref name="numColsRight"/>.
        /// Input checking is the responsibility of the caller.
        /// </summary>
        /// <param name="numRowsLeft">The number of rows of <paramref name="matrixLeft"/>. Constraint:
        ///     <paramref name="numRowsLeft"/> == <paramref name="numRowsRight"/>.</param>
        /// <param name="numColsLeft">The number of columns of <paramref name="matrixLeft"/>.</param>
        /// <param name="matrixLeft">The entries of the left matrix in column major layout. It will not be modified.</param>
        /// <param name="numRowsRight">The number of rows of <paramref name="matrixRight"/>. Constraint:
        ///     <paramref name="numRowsLeft"/> == <paramref name="numRowsRight"/>.</param>
        /// <param name="numColsRight">The number of columns of <paramref name="matrixRight"/>.</param>
        /// <param name="matrixRight">The entries of the right matrix in column major layout. It will not be modified.</param>
        /// <returns></returns>
        internal static double[] JoinHorizontally(int numRowsLeft, int numColsLeft, double[] matrixLeft,
            int numRowsRight, int numColsRight, double[] matrixRight)
        {
            double[] result = new double[matrixLeft.Length + matrixRight.Length];
            // In result all columns of matrixLeft will be contiguous and followed by all columns of matrixRight, also contiguous
            Array.Copy(matrixLeft, result, matrixLeft.Length); 
            Array.Copy(matrixRight, 0, result, matrixLeft.Length, matrixRight.Length);
            return result;
        }

        /// <summary>
        /// result = [<paramref name="matrixUp"/>; <paramref name="matrixDown"/>] (Matlab notation). Dimensions of result:
        /// (<paramref name="numRowsUp"/> + <paramref name="numRowsDown"/>) -by-
        /// (<paramref name="numColsUp"/> = <paramref name="numColsDown"/>.
        /// Input checking is the responsibility of the caller.
        /// </summary>
        /// <param name="numRowsUp">The number of rows of <paramref name="matrixUp"/>.</param>
        /// <param name="numColsUp">The number of columns of <paramref name="matrixUp"/>. Constraint:
        ///     <paramref name="numColsUp"/> == <paramref name="numColsDown"/>.</param>
        /// <param name="matrixUp">The entries of the top matrix in column major layout. It will not be modified.</param>
        /// <param name="numRowsDown">The number of rows of <paramref name="matrixDown"/>.</param>
        /// <param name="numColsDown">The number of columns of <paramref name="matrixDown"/>. Constraint:
        ///     <paramref name="numColsUp"/> == <paramref name="numColsDown"/>.</param>
        /// <param name="matrixDown">The entries of the bottom matrix in column major layout. It will not be modified.</param>
        /// <returns></returns>
        internal static double[] JoinVertically(int numRowsUp, int numColsUp, double[] matrixUp,
            int numRowsDown, int numColsDown, double[] matrixDown)
        {
            double[] result = new double[matrixUp.Length + matrixDown.Length];
            int numRowsResult = numRowsUp + numRowsDown;
            // Copy one contiguous column at a time. 
            for (int j = 0; j < numColsUp; ++j) // TODO: processing matrixUp and then matrixDown might be faster
            {
                Array.Copy(matrixUp, j * numRowsUp, result, j * numRowsResult, numRowsUp);
                Array.Copy(matrixDown, j * numRowsDown, result, j * numRowsResult + numRowsUp, numRowsDown);
            }
            return result;
        }

        /// <summary>
        /// Add zero rows to the bottom of the matrix.
        /// </summary>
        /// <param name="numRowsOld">The number of rows of <paramref name="matrix"/>.</param>
        /// <param name="numCols">The number of columns of <paramref name="matrix"/></param>
        /// <param name="numRowsNew">The number of columns after the operation. It must be 
        ///     <paramref name="numRowsNew"/> &gt; <paramref name="numRowsOld"/>.</param>
        /// <param name="matrix">The matrix in full column major layout.</param>
        /// <returns></returns>
        internal static double[] IncreaseRows(int numRowsOld, int numCols, int numRowsNew, double[] matrix)
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

        /// <summary>
        /// Exchanges the rows and columns of a square matrix, according to the provided permutation.
        /// Input checking is the responsibility of the caller.
        /// </summary>
        /// <param name="order">The number of rows = number of columns of <paramref name="matrix"/>.</param>
        /// <param name="matrix">The matrix in full column major layout.</param>
        /// <param name="permutation">A map between the old and new indices: 
        ///     newIndex = <paramref name="permutation"/>[oldIndex].</param>
        /// <returns></returns>
        public static double[] ReorderOldToNew(int order, double[] matrix, IReadOnlyList<int> permutation)
        {
            double[] reordered = new double[matrix.Length];
            for (int j = 0; j < order; ++j)
            {
                int newCol = permutation[j];
                for (int i = 0; i < order; ++i) reordered[newCol * order + permutation[i]] = matrix[j * order + i];
            }
            return reordered;
        }

        /// <summary>
        /// Exchanges the rows and columns of a square matrix, according to the provided permutation.
        /// Input checking is the responsibility of the caller.
        /// </summary>
        /// <param name="order">The number of rows = number of columns of <paramref name="matrix"/>.</param>
        /// <param name="matrix">The matrix in full column major layout.</param>
        /// <param name="permutation">A map between the old and new indices: 
        ///     oldIndex = <paramref name="permutation"/>[newIndex].</param>
        /// <returns></returns>
        public static double[] ReorderNewToOld(int order, double[] matrix, IReadOnlyList<int> permutation)
        {
            double[] reordered = new double[matrix.Length];
            for (int j = 0; j < order; ++j)
            {
                int oldCol = permutation[j];
                for (int i = 0; i < order; ++i) reordered[j * order + i] = matrix[oldCol * order + permutation[i]];
            }
            return reordered;
        }

        /// <summary>
        /// Add zero columns to the right of the matrix or remove existing ones.
        /// </summary>
        /// <param name="numRows">The number of rows in <paramref name="matrix"/>.</param>
        /// <param name="numColsOld">The number of columns in <paramref name="matrix"/>.</param>
        /// <param name="numColsNew">The number of columns after the operation.</param>
        /// <param name="matrix">The matrix in full column major layout.</param>
        /// <returns></returns>
        internal static double[] ResizeCols(int numRows, int numColsOld, int numColsNew, double[] matrix)
        {
            // We need to copy the common columns. Since the matrix is column major, they are contiguous
            double[] result = new double[numRows * numColsNew];
            int numColsMin = (numColsNew > numColsOld) ? numColsOld : numColsNew;
            Array.Copy(matrix, result, matrix.Length);
            return result;
        }

        /// <summary>
        /// Sets a column <paramref name="colIdx"/> of <paramref name="matrix"/> to <paramref name="newCol"/>.
        /// Input checking is the responsibility of the caller.
        /// </summary>
        /// <param name="numRows">The number of matrix rows.</param>
        /// <param name="numCols">The number of matrix rows columns.</param>
        /// <param name="matrix">The matrix entries in column major layout. This array will be modified.</param>
        /// <param name="colIdx">The index of the column to alter: 
        ///     0 &lt;= <paramref name="colIdx"/> &lt; <paramref name="numCols"/></param>
        /// <param name="rowStart">The (inclusive) index of the row of column <paramref name="colIdx"/> after which 
        ///     <paramref name="newCol"/> will be copied.</param>
        /// <param name="newCol">The new entries of column <paramref name="colIdx"/>. Its length must be &lt;= 
        ///     <paramref name="numRows"/></param>
        internal static void SetCol(int numRows, int numCols, double[] matrix, int colIdx, int rowStart, double[] newCol)
        {
            // Column entries are contiguous:
            Array.Copy(newCol, 0, matrix, colIdx * numRows + rowStart, newCol.Length);
        }

        /// <summary>
        /// Sets a row <paramref name="rowIdx"/> of <paramref name="matrix"/> to <paramref name="newRow"/>.
        /// Input checking is the responsibility of the caller.
        /// </summary>
        /// <param name="numRows">The number of matrix rows.</param>
        /// <param name="numCols">The number of matrix rows columns.</param>
        /// <param name="matrix">The matrix entries in column major layout. This array will be modified.</param>
        /// <param name="rowIdx">The index of the row to alter: 
        ///     0 &lt;= <paramref name="rowIdx"/> &lt; <paramref name="numRows"/></param>
        /// <param name="colStart">The (inclusive) index of the column of row <paramref name="rowIdx"/> after which 
        ///     <paramref name="newRow"/> will be copied.</param>
        /// <param name="newRow">The new entries of row <paramref name="rowIdx"/>. Its length must be &lt;=  
        ///     <paramref name="numCols"/></param>
        internal static void SetRow(int numRows, int numCols, double[] matrix, int rowIdx, int colStart, double[] newRow)
        {
            // Row entries are not contiguous
            for (int j = colStart; j < colStart + newRow.Length; ++j) matrix[j * numRows + rowIdx] = newRow[j];
        }

        /// <summary>
        /// Sets the entries of <paramref name="matrixBig"/> between rows 
        /// [<paramref name="rowStartBig"/>, <paramref name="rowStartBig"/> + <paramref name="numRowsSub"/>) 
        /// and between columns 
        /// [<paramref name="colStartBig"/>, <paramref name="colStartBig"/> + <paramref name="numColsSub"/>) to the entries of
        /// <paramref name="matrixSub"/>.
        /// Input checking is the responsibility of the caller.
        /// </summary>
        /// <param name="numRowsBig">The number of rows of the big matrix <paramref name="matrixBig"/>.</param>
        /// <param name="numColsBig">The number of columns of the big matrix <paramref name="matrixBig"/>.</param>
        /// <param name="matrixBig">The entries of the big matrix in column major layout. This array will be modified.</param>
        /// <param name="rowStartBig">The index of the first row of <paramref name="matrixBig"/> that will be altered. Constraint:
        ///     <paramref name="rowStartBig"/> + <paramref name="numRowsSub"/> &lt;= <paramref name="numRowsBig"/>.</param>
        /// <param name="colStartBig">The index of the first column of <paramref name="matrixBig"/> that will be altered. Constraint:
        ///     <paramref name="colStartBig"/> + <paramref name="numColsSub"/> &lt;= <paramref name="numColsBig"/>.</param>
        /// <param name="numRowsSub">The number of rows of the submatrix <paramref name="matrixSub"/>. Constraint:
        ///     <paramref name="rowStartBig"/> + <paramref name="numRowsSub"/> &lt;= <paramref name="numRowsBig"/>.</param>
        /// <param name="numColsSub">The number of columns of the big matrix <paramref name="matrixSub"/>. Constraint:
        ///     <paramref name="colStartBig"/> + <paramref name="numColsSub"/> &lt;= <paramref name="numColsBig"/>.</param>
        /// <param name="matrixSub">The entries of the submatrix in column major layout. This array will not be modified.</param>
        internal static void SetSubmatrix(int numRowsBig, int numColsBig, double[] matrixBig, int rowStartBig, int colStartBig,
            int numRowsSub, int numColsSub, double[] matrixSub)
        {
            // Optimization case: if the submatrix has exactly as many rows as the big matrix and there is now row offset, then 
            // the entries of the big matrix that will be set are contiguous.
            if ((numRowsBig == numRowsSub) && (rowStartBig == 0))
            {
                Array.Copy(matrixSub, 0, matrixBig, colStartBig * numRowsBig, numRowsSub * numColsSub);
            }
            else
            {
                for (int j = 0; j < numColsSub; ++j)
                {
                    // The numRowsSub entries of column j are contigous.
                    Array.Copy(matrixSub, j * numRowsSub, matrixBig, (colStartBig + j) * numRowsBig + rowStartBig, numRowsSub);
                    // TODO: if the rows of submatrix are very few, the overhead of Array.Copy() may be worse than directly 
                    // setting the entries.
                }
            }
        }
    }
}
