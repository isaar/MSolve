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

        public static readonly double[] x10 = { 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9, 10.1 };
        public static readonly double[] b5_1 = { 5.5, 0.0, -2.2, 0.0, -8.8 };
        public static readonly double[] b5_2 = { -1.1, -1.1, -1.1, -1.1, -1.1 };
        public static readonly double[] x5 = { 1.1, 2.2, 3.3, 4.4, 5.5 };
        public static readonly double[] b10_1 = { 0.0,-3.3, 0.0, 0.0, 1.1, 0.0, 0.0,-5.5, 0.0, 0.0 };
        public static readonly double[] b10_2 = { 1.1, 1.1, 1.1, 1.1, 1.1, -5.5, 0.0, 0.0, 0.0, 0.0, };

        public static void CheckIndexer()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            SignedBooleanMatrix B1 = CreateMatrix1();
            SignedBooleanMatrix B2 = CreateMatrix2();
            comparer.CheckMatrixEquality(matrix1, B1.CopyToArray2D());
            comparer.CheckMatrixEquality(matrix2, B2.CopyToArray2D());
        }

        public static void CheckMatrixVectorMultiplication()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            SignedBooleanMatrix B1 = CreateMatrix1();
            SignedBooleanMatrix B2 = CreateMatrix2();
            Vector B1timesX1 = B1.MultiplyRight(Vector.CreateFromArray(x10, true), false);
            Vector B2timesX1 = B2.MultiplyRight(Vector.CreateFromArray(x10, true), false);
            Vector transB1timesX2 = B1.MultiplyRight(Vector.CreateFromArray(x5, true), true);
            Vector transB2timesX2 = B2.MultiplyRight(Vector.CreateFromArray(x5, true), true);
            comparer.CheckMatrixVectorMult(matrix1, x10, b5_1, B1timesX1.CopyToArray());
            comparer.CheckMatrixVectorMult(matrix2, x10, b5_2, B2timesX1.CopyToArray());
            comparer.CheckMatrixVectorMult(matrix1, x5, b10_1, transB1timesX2.CopyToArray());
            comparer.CheckMatrixVectorMult(matrix2, x5, b10_2, transB2timesX2.CopyToArray());
        }

        public static void WriteToConsole()
        {
            SignedBooleanMatrix B1 = CreateMatrix1();
            SignedBooleanMatrix B2 = CreateMatrix2();
            B1.WriteToConsole();
            B2.WriteToConsole();
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
