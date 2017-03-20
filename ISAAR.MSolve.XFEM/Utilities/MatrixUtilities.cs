using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.XFEM.Utilities
{
    /// <summary>
    /// TODO: Some of these could be extension methods, but it effort should be spent into implementing them in the  
    /// Matrix classes directly. If they indeed only need IMatrix, then they should be implemented as extensions.
    /// </summary>
    static class MatrixUtilities
    {
        public static void AddPartialToTotalMatrix(Matrix2D<double> partialMatrix, Matrix2D<double> totalMatrix)
        {
            Debug.Assert(partialMatrix.Rows == totalMatrix.Rows, "Non matching rows");
            Debug.Assert(partialMatrix.Columns == totalMatrix.Columns, "Non matching columns");
            for (int row = 0; row < totalMatrix.Rows; ++row)
            {
                for (int col = 0; col < totalMatrix.Columns; ++col)
                {
                    totalMatrix[row, col] += partialMatrix[row, col];
                }
            }
        }

        public static void AddPartialToSymmetricTotalMatrix(Matrix2D<double> partialMatrix,
            SymmetricMatrix2D<double> totalMatrix)
        {
            Debug.Assert(partialMatrix.Rows == totalMatrix.Rows, "Non matching rows");
            Debug.Assert(partialMatrix.Columns == totalMatrix.Columns, "Non matching columns");
            for (int row = 0; row < totalMatrix.Rows; ++row)
            {
                for (int col = row; col < totalMatrix.Columns; ++col)
                {
                    totalMatrix[row, col] += partialMatrix[row, col];
                }
            }
        }

        public static void PrintDense(IMatrix2D<double> matrix)
        {
            StringBuilder builder = new StringBuilder();
            for (int row = 0; row < matrix.Rows; ++row)
            {
                for (int col = 0; col < matrix.Columns; ++col)
                {
                    builder.Append(matrix[row, col]);
                    builder.Append(' ');
                }
                builder.Append("\n");
            }
            builder.Append("\n");
            Console.Write(builder.ToString());
        }

        public static void PrintUpperTriangleDense(IMatrix2D<double> matrix)
        {
            StringBuilder builder = new StringBuilder();
            for (int row = 0; row < matrix.Rows; ++row)
            {
                for (int col = 0; col < matrix.Columns; ++col)
                {
                    if (col < row) builder.Append(0.0);
                    else builder.Append(matrix[row, col]);
                    builder.Append(' ');
                }
                builder.Append("\n");
            }
            builder.Append("\n");
            Console.Write(builder.ToString());
        }

        /// <summary>
        /// The output is a dense matrix
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="oldToNewIndices"></param>
        /// <returns></returns>
        public static Matrix2D<double> Reorder(IMatrix2D<double> matrix, int[] oldToNewIndices)
        {
            int order = matrix.Rows;
            if (order != matrix.Columns) throw new ArgumentException("The matrix must be square");
            if (order != oldToNewIndices.Length)
                throw new ArgumentException("Mismatch in the dimensions of the matrix and the permutation array");

            var newMatrix = new Matrix2D<double>(order, order);
            for (int row = 0; row < order; ++row)
            {
                int newRow = oldToNewIndices[row];
                for (int col = 0; col < order; ++col)
                {
                    newMatrix[newRow, oldToNewIndices[col]] = matrix[row, col];
                }
            }
            return newMatrix;
        }
    }
}
