using System;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    static class Preconditions
    {
        public static void CheckIndices(IIndexable2D matrix, int rowIdx, int colIdx)
        {
            if ((rowIdx < 0) || (rowIdx >= matrix.NumRows))
            {
                throw new IndexOutOfRangeException($"Row index {rowIdx} was outside the range"
                    + $" [0, number_of_rows={matrix.NumRows})");
            }
            if ((colIdx < 0) || (colIdx >= matrix.NumColumns))
            {
                throw new IndexOutOfRangeException($"Column index {colIdx} was outside the range"
                    + $" [0, number_of_columns={matrix.NumColumns})");
            }
        }

        public static void CheckVectorDimensions(IVectorView vector1, IVectorView vector2)
        {
            if (vector1.Length != vector2.Length)
            {
                string message = string.Format("Vector1 has length of {0}, while vector2 has length of {1}", 
                    vector1.Length, vector2.Length);
                throw new NonMatchingDimensionsException(message);
            }
               
        }

        public static bool AreSameMatrixDimensions(IIndexable2D matrix1, IIndexable2D matrix2)
        {
            if ((matrix1.NumRows != matrix2.NumRows) || (matrix1.NumColumns != matrix2.NumColumns)) return false;
            else return true;
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
    }
}
