using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    static class Preconditions
    {
        public static void CheckVectorDimensions(IVectorView vector1, IVectorView vector2)
        {
            if (vector1.Length != vector2.Length)
            {
                string message = string.Format("Vector1 has length of {0}, while vector2 has length of {1}", 
                    vector1.Length, vector2.Length);
                throw new NonMatchingDimensionsException(message);
            }
               
        }

        public static void CheckSameMatrixDimensions(IMatrixView matrix1, IMatrixView matrix2)
        {
            if ((matrix1.NumRows != matrix2.NumRows) != (matrix1.NumColumns != matrix2.NumColumns))
            {
                string message = string.Format(
                    "Matrix1 has dimensions ({0}x{1}), while matrix2 has dimensions ({2}x{3})",
                    matrix1.NumRows, matrix1.NumColumns, matrix2.NumRows, matrix2.NumColumns);
                throw new NonMatchingDimensionsException(message);
            }
        }

        public static void CheckMultiplicationDimensions(IMatrixView leftMatrix, IMatrixView rightMatrix)
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

        public static void CheckSystemSolutionDimensions(IMatrixView matrix, IVectorView rhsVector)
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
