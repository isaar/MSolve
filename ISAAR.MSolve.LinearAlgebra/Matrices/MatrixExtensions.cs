using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
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

        public static IMatrixView Scale(this IMatrixView matrix, double scalar)
        {
            return matrix.DoToAllEntries(x => scalar * x);
        }

        public static void ScaleIntoThis(this IMatrix matrix, double scalar)
        {
            matrix.DoToAllEntriesIntoThis(x => scalar * x);
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
    }
}
