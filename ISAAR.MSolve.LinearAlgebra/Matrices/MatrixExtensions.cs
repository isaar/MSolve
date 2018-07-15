using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// If one of the following methods is also implemented in a concrete class, use that one (if possible), as it will be far 
    /// more efficient.
    /// </summary>
    public static class MatrixExtensions
    {
        /// <summary>
        /// result[i, j] = matrix1[i, j] + matrix2[i, j]
        /// </summary>
        /// <param name="matrix1"></param>
        /// <param name="matrix2"></param>
        /// <returns></returns>
        public static IMatrixView Add(this IMatrixView matrix1, IMatrixView matrix2)
        {
            return matrix1.Axpy(matrix2, 1.0);
        }

        /// <summary>
        /// matrix1 = matrix1 + matrix2
        /// </summary>
        /// <param name="matrix1"></param>
        /// <param name="matrix2"></param>
        public static void AddIntoThis(this IMatrix matrix1, IMatrixView matrix2)
        {
            matrix1.AxpyIntoThis(matrix2, 1.0);
        }

        public static Vector GetDiagonal(this IMatrixView matrix)
        {
            return Vector.CreateFromArray(matrix.GetDiagonalAsArray(), false);
        }

        public static double[] GetDiagonalAsArray(this IMatrixView matrix)
        {
            Preconditions.CheckSquare(matrix);
            double[] diag = new double[matrix.NumRows];
            for (int i = 0; i < matrix.NumRows; ++i) diag[i] = matrix[i, i];
            return diag;
        }

        public static bool IsSymmetric(this IIndexable2D matrix, double tolerance = double.Epsilon)
        {
            var comparer = new ValueComparer(tolerance);
            if (matrix.NumRows != matrix.NumColumns) return false;
            for (int i = 0; i < matrix.NumRows; ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    if (!comparer.AreEqual(matrix[i, j], matrix[j, i])) return false;
                }
            }
            return true;
        }

        public static Matrix Reorder(this IIndexable2D matrix, IReadOnlyList<int> permutation, bool oldToNew)
        {
            int order = matrix.NumRows;
            if (matrix.NumColumns != order)
            {
                throw new NonMatchingDimensionsException("This operation works on square matrices only.");
            }
            if (permutation.Count != order)
            {
                throw new NonMatchingDimensionsException($"This matrix has order = {order}, while the permutation vector"
                    + $" has {permutation.Count} entries.");
            }
            var reordered = Matrix.CreateZero(order, order);
            if (oldToNew)
            {
                for (int j = 0; j < order; ++j)
                {
                    int newCol = permutation[j];
                    for (int i = 0; i < order; ++i) reordered[permutation[i], newCol] = matrix[i, j];
                }
            }
            else
            {
                for (int j = 0; j < order; ++j)
                {
                    int oldCol = permutation[j];
                    for (int i = 0; i < order; ++i) reordered[i, j] = matrix[permutation[i], oldCol];
                }
            }
            return reordered;
        }

        /// <summary>
        /// result[i, j] = matrix1[i, j] - matrix2[i, j]
        /// </summary>
        /// <param name="matrix1"></param>
        /// <param name="matrix2"></param>
        /// <returns></returns>
        public static IMatrixView Subtract(this IMatrixView matrix1, IMatrixView matrix2)
        {
            return matrix1.Axpy(matrix2, - 1.0);
        }

        /// <summary>
        /// matrix1 = matrix1 - matrix2
        /// </summary>
        /// <param name="matrix1"></param>
        /// <param name="matrix2"></param>
        public static void SubtractIntoThis(this IMatrix matrix1, IMatrixView matrix2)
        {
            matrix1.AxpyIntoThis(matrix2 ,- 1.0);
        }

        public static Numerical.LinearAlgebra.Matrix2D MultiplyTransposeThisTimesOtherTimesThis(this CSCMatrix csc, 
            Numerical.LinearAlgebra.Matrix2D other)
        {
            var otherCopy = Matrix.CreateFromLegacyMatrix(other);
            Matrix otherTimesThis = csc.MultiplyLeft(otherCopy, false, false);
            Matrix result = csc.MultiplyRight(otherTimesThis, true, false);
            return new Numerical.LinearAlgebra.Matrix2D(result.CopyToArray2D());
        }

        public static double[] MultiplyRight(this CSCMatrix csc, double[] vector, bool tranposeThis)
        {
            var asVector = Vector.CreateFromArray(vector, false);
            return csc.MultiplyRight(asVector, tranposeThis).InternalData;
        }
    }
}
