using System;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Commons
{
    /// <summary>
    /// Checks for the compatibility of vector or matrix dimensions for various operations and for the indices. Also templates 
    /// for exception messages.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class Preconditions
    {
        public static bool AreSameMatrixDimensions(IIndexable2D matrix1, IIndexable2D matrix2)
        {
            if ((matrix1.NumRows != matrix2.NumRows) || (matrix1.NumColumns != matrix2.NumColumns)) return false;
            else return true;
        }

        public static void CheckIndex1D(IVectorView vector, int idx)
        {
            if ((idx < 0) || (idx >= vector.Length))
            {
                throw new IndexOutOfRangeException($"Cannot access index {idx} in a vector with length = {vector.Length}");
            }
        }

        public static void CheckIndexCol(IIndexable2D matrix, int colIdx)
        {
            if ((colIdx < 0) || (colIdx >= matrix.NumColumns))
            {
                throw new IndexOutOfRangeException($"Cannot access column {colIdx} in a"
                    + $" {matrix.NumRows}-by-{matrix.NumColumns} matrix");
            }
        }

        public static void CheckIndexRow(IIndexable2D matrix, int rowIdx)
        {
            if ((rowIdx < 0) || (rowIdx >= matrix.NumRows))
            {
                throw new IndexOutOfRangeException($"Cannot access row {rowIdx} in a"
                    + $" {matrix.NumRows}-by-{matrix.NumColumns} matrix");
            }
        }

        public static void CheckIndices(IIndexable2D matrix, int rowIdx, int colIdx)
        {
            if ((rowIdx < 0) || (rowIdx >= matrix.NumRows))
            {
                throw new IndexOutOfRangeException($"Cannot access row {rowIdx} in a"
                    + $" {matrix.NumRows}-by-{matrix.NumColumns} matrix");
            }
            if ((colIdx < 0) || (colIdx >= matrix.NumColumns))
            {
                throw new IndexOutOfRangeException($"Cannot access column {colIdx} in a"
                    + $" {matrix.NumRows}-by-{matrix.NumColumns} matrix");
            }
        }

        public static void CheckMultiplicationDimensions(IIndexable2D leftMatrix, IIndexable2D rightMatrix)
        {
            if (leftMatrix.NumColumns != rightMatrix.NumRows)
            {
                string message = string.Format(
                    "Left matrix has dimensions ({0}x{1}), while right matrix has dimensions ({2}x{3})",
                    leftMatrix.NumRows, leftMatrix.NumColumns, rightMatrix.NumRows, rightMatrix.NumColumns);
                throw new NonMatchingDimensionsException(message);
            }
        }

        //This might be better for the methods that support transposition of the left or right matrix
        public static void CheckMultiplicationDimensions(int leftColumns, int rightRows)
        {
            if (leftColumns != rightRows)
            {
                string message = $"Left matrix/vector has {leftColumns} columns, while right matrix/vector has {rightRows} rows";
                throw new NonMatchingDimensionsException(message);
            }
        }

        public static void CheckMultiplicationDimensionsSection(IIndexable2D matrixLeft, IVectorView vectorRight,
            int vectorStart, IVectorView result, int resultStart)
        {
            if (vectorStart + matrixLeft.NumColumns > vectorRight.Length) throw new NonMatchingDimensionsException(
                $"The multiplied vector's length = {vectorRight.Length} must be at least as large as the start index =" +
                $" {vectorStart} + the matrix' columns = {matrixLeft.NumColumns}");
            if (vectorStart + matrixLeft.NumRows > result.Length) throw new NonMatchingDimensionsException(
                $"The result vector's length = {result.Length} must be at least as large as the start index =" +
                $" {resultStart} + the matrix' rows = {matrixLeft.NumRows}");
        }

        public static void CheckSameColDimension(IIndexable2D matrix, IVectorView vector)
        {
            if (matrix.NumColumns != vector.Length)
            {
                string message = string.Format(
                    "Matrix has dimensions ({0}x{1}), while vector has dimensions ({2}x{1})",
                    matrix.NumRows, matrix.NumColumns, vector.Length);
                throw new NonMatchingDimensionsException(message);
            }
        }

        public static void CheckSameColDimension(IIndexable2D matrix1, IIndexable2D matrix2)
        {
            if (matrix1.NumColumns != matrix2.NumColumns)
            {
                string message = string.Format(
                    "Matrix1 has dimensions ({0}x{1}), while matrix2 has dimensions ({2}x{3})",
                    matrix1.NumRows, matrix1.NumColumns, matrix2.NumRows, matrix2.NumColumns);
                throw new NonMatchingDimensionsException(message);
            }
        }

        public static void CheckSameMatrixDimensions(IIndexable2D matrix1, IIndexable2D matrix2)
        {
            if ((matrix1.NumRows != matrix2.NumRows) || (matrix1.NumColumns != matrix2.NumColumns))
            {
                string message = string.Format(
                    "Matrix1 has dimensions ({0}x{1}), while matrix2 has dimensions ({2}x{3})",
                    matrix1.NumRows, matrix1.NumColumns, matrix2.NumRows, matrix2.NumColumns);
                throw new NonMatchingDimensionsException(message);
            }
        }

        public static void CheckSameRowDimension(IIndexable2D matrix, IVectorView vector)
        {
            if (matrix.NumRows != vector.Length)
            {
                string message = string.Format(
                    "Matrix has dimensions ({0}x{1}), while vector has dimensions ({2}x{1})",
                    matrix.NumRows, matrix.NumColumns, vector.Length);
                throw new NonMatchingDimensionsException(message);
            }
        }

        public static void CheckSameRowDimension(IIndexable2D matrix1, IIndexable2D matrix2)
        {
            if (matrix1.NumRows != matrix2.NumRows)
            {
                string message = string.Format(
                    "Matrix1 has dimensions ({0}x{1}), while matrix2 has dimensions ({2}x{3})",
                    matrix1.NumRows, matrix1.NumColumns, matrix2.NumRows, matrix2.NumColumns);
                throw new NonMatchingDimensionsException(message);
            }
        }

        public static void CheckSquare(IIndexable2D matrix)
        {
            if (matrix.NumRows != matrix.NumColumns) throw new NonMatchingDimensionsException(
                $"The matrix must be square, but was {matrix.NumRows}-by-{matrix.NumColumns}");
        }

        public static void CheckSystemSolutionDimensions(IIndexable2D matrix, IVectorView rhsVector)
        {
            if (matrix.NumRows != rhsVector.Length)
            {
                string message = string.Format(
                    "Matrix has dimensions ({0}x{1}), while the right hand side vector has dimensions ({2}x1)",
                    matrix.NumRows, matrix.NumColumns, rhsVector.Length);
                throw new NonMatchingDimensionsException(message);
            }
        }

        public static void CheckSystemSolutionDimensions(int matrixRows, int matrixColumns, int rhsVectorLength)
        {
            if (matrixRows != rhsVectorLength)
            {
                string message = $"Matrix has dimensions ({matrixRows}x{matrixColumns}), while the right hand side vector has"
                    + $" dimensions ({rhsVectorLength}x1)";
                throw new NonMatchingDimensionsException(message);
            }
        }

        public static void CheckVectorDimensions(IIndexable1D vector1, IIndexable1D vector2)
        {
            if (vector1.Length != vector2.Length)
            {
                string message = string.Format("Vector1 has length of {0}, while vector2 has length of {1}",
                    vector1.Length, vector2.Length);
                throw new NonMatchingDimensionsException(message);
            }

        }
    }
}
