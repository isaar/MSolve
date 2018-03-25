using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    public static class MatrixExtensions
    {
        /// <summary>
        /// result = matrix1 + matrix2
        /// </summary>
        /// <param name="matrix1"></param>
        /// <param name="matrix2"></param>
        /// <returns></returns>
        public static Matrix Add(this Matrix matrix1, Matrix matrix2)
        {
            return matrix1.Axpy(1.0, matrix2);
        }

        /// <summary>
        /// matrix1 = matrix1 + matrix2
        /// </summary>
        /// <param name="matrix1"></param>
        /// <param name="matrix2"></param>
        public static void AddIntoThis(this Matrix matrix1, Matrix matrix2)
        {
            matrix1.AxpyIntoThis(1.0, matrix2);
        }

        /// <summary>
        /// result = matrix1 - matrix2
        /// </summary>
        /// <param name="matrix1"></param>
        /// <param name="matrix2"></param>
        /// <returns></returns>
        public static Matrix Subtract(this Matrix matrix1, Matrix matrix2)
        {
            return matrix1.Axpy(-1.0, matrix2);
        }

        /// <summary>
        /// matrix1 = matrix1 - matrix2
        /// </summary>
        /// <param name="matrix1"></param>
        /// <param name="matrix2"></param>
        public static void SubtractIntoThis(this Matrix matrix1, Matrix matrix2)
        {
            matrix1.AxpyIntoThis(-1.0, matrix2);
        }

        public static void WriteToConsole(this IIndexable2D matrix, Array2DFormatting format = null)
        {
            if (format == null) format = Array2DFormatting.Brackets;
            Console.Write(format.ArrayStart);

            // First row
            Console.Write(format.RowSeparator + format.RowStart);
            Console.Write(matrix[0, 0]);
            for (int j = 1; j < matrix.NumColumns; ++j)
            {
                Console.Write(format.ColSeparator + matrix[0, j]);
            }
            Console.Write(format.RowEnd);

            // Subsequent rows
            for (int i = 1; i < matrix.NumRows; ++i)
            {
                Console.Write(format.RowSeparator + format.RowStart);
                Console.Write(matrix[i, 0]);
                for (int j = 1; j < matrix.NumColumns; ++j)
                {
                    Console.Write(format.ColSeparator + matrix[i, j]);
                }
                Console.Write(format.RowEnd);
            }
            Console.Write(format.RowSeparator + format.ArrayEnd);
        }
    }
}
