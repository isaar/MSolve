using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public static class MatrixExtensions
    {
        public static void PrintMatrix2D(this IMatrix2D matrix, string filePath)
        {
            using (var writer = new StreamWriter(filePath))
            {
                for (int i = 0; i < matrix.Rows; ++i)
                {
                    for (int j = 0; j < matrix.Columns; ++j)
                    {
                        writer.Write(matrix[i, j]);
                        writer.Write(" ");
                    }
                    writer.WriteLine();
                }
                writer.Flush();
            }
        }
    }
}
