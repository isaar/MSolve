using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    internal static class MultiDiagonalMatrices
    {
        internal static double[,] CreateSymmetric(int matrixOrder, int[] diagonalOrdinals)
        {
            var matrix = new double[matrixOrder, matrixOrder];
            diagonalOrdinals = diagonalOrdinals.Where(d => d != 0).ToArray(); // Ignore the 0 diagonal since it is the main one.

            double minValueMain = matrixOrder;
            for (int i = 0; i < matrixOrder; ++i) matrix[i, i] = minValueMain + i;

            foreach (int d in diagonalOrdinals)
            {
                // Divide over a function of the diagonal and the number of diagonals to keep the matrix pos def.
                double minValue = minValueMain / (diagonalOrdinals.Length + d); 
                for (int i = 0; i < matrixOrder - d; ++i)
                {
                    matrix[i, i + d] = minValue + i;
                    matrix[i + d, i] = minValue + i;
                }
            }

            //var writer = new LinearAlgebra.Output.Array2DWriter();
            //writer.WriteToFile(matrix, @"C:\Users\Serafeim\Desktop\output.txt");

            return matrix;
        }
    }
}
