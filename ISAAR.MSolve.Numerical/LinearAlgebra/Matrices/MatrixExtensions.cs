using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
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
