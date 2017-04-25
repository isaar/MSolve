using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.LinearAlgebra
{
    /// <summary>
    /// TODO: Some of these could be extension methods, but effort should be spent into implementing them in the  
    /// Matrix classes directly. If they indeed only need IMatrix, then they should be implemented as extensions.
    /// </summary>
    static class MatrixUtilities
    {
        public static void AddPartialToTotalMatrix(Matrix2D partialMatrix, Matrix2D totalMatrix)
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

        public static void AddPartialToSymmetricTotalMatrix(Matrix2D partialMatrix,
            SymmetricMatrix2D totalMatrix)
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

        public static void PrintDense(IMatrix2D matrix)
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

        public static void PrintUpperTriangleDense(IMatrix2D matrix)
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
        public static Matrix2D Reorder(IMatrix2D matrix, int[] oldToNewIndices)
        {
            int order = matrix.Rows;
            if (order != matrix.Columns) throw new ArgumentException("The matrix must be square");
            if (order != oldToNewIndices.Length) throw new NonMatchingDimensionsException(
                "Mismatch in the dimensions of the matrix and the permutation array");

            var newMatrix = new Matrix2D(order, order);
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

        public static double[] Substract(double[] vector1, double[] vector2)
        {
            if (vector1.Length != vector2.Length) throw new NonMatchingDimensionsException(
                "Vector 1 length  = " + vector1.Length + ", but vector 2 length = " + vector2.Length);

            double[] result = new double[vector1.Length];
            for (int i = 0; i < vector1.Length; ++i) result[i] = vector1[i] - vector2[i];
            return result;
        }
    }
}
