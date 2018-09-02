using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices
{
    class SignedBoolean
    {
        public const int numRows = 5;
        public const int numCols = 10;

        public static readonly double[,] matrix1 =
        {
            { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
            { 0,-1, 0, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 0, 0,-1, 0, 0 }
        };

        public static readonly double[,] matrix2 =
        {
            { 1,-1, 0, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 1,-1, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 1,-1, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 1,-1, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 1,-1, 0, 0, 0, 0 }
        };

        public static void WriteToConsole()
        {
            var writer = new BooleanMatrixWriter(true);
            SignedBooleanMatrix B1 = CreateMatrix1();
            SignedBooleanMatrix B2 = CreateMatrix2();
            writer.WriteToConsole(B1);
            writer.WriteToConsole(B2);
        }

        private static SignedBooleanMatrix CreateMatrix1()
        {
            var B1 = new SignedBooleanMatrix(numRows, numCols);
            for (int i = 0; i < numRows; ++i)
            {
                for (int j = 0; j < numCols; ++j)
                {
                    if (matrix1[i, j] == 1.0) B1.AddEntry(i, j, true);
                    else if (matrix1[i, j] == -1.0) B1.AddEntry(i, j, false);
                }
            }
            return B1;
        }

        private static SignedBooleanMatrix CreateMatrix2()
        {
            var B2 = new SignedBooleanMatrix(numRows, numCols);
            for (int i = 0; i < numRows; ++i)
            {
                for (int j = 0; j < numCols; ++j)
                {
                    if (matrix2[i, j] == 1.0) B2.AddEntry(i, j, true);
                    else if (matrix2[i, j] == -1.0) B2.AddEntry(i, j, false);
                }
            }
            return B2;
        }
    }
}
