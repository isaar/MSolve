using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Reordering;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices
{
    class MatrixAssembly
    {
        public const int subOrder = 4;
        public const int globalOrder = 8;

        public static readonly double[,] subMatrix1 =
        {
            { 20.1,  1.1,  2.1,  3.1 }, 
            {  1.1, 20.2,  2.2,  3.2 },
            {  2.1,  2.2, 20.3,  3.3 },
            {  3.1,  3.2,  3.3, 20.4 }
        };

        public static readonly double[,] subMatrix2 =
        {
            { 30.1,  1.1,  2.1,  3.1 },
            {  1.1, 30.2,  2.2,  3.2 },
            {  2.1,  2.2, 30.3,  3.3 },
            {  3.1,  3.2,  3.3, 30.4 }
        };

        public static readonly double[,] subMatrix3 =
        {
            { 40.1,  1.1,  2.1,  3.1 },
            {  1.1, 40.2,  2.2,  3.2 },
            {  2.1,  2.2, 40.3,  3.3 },
            {  3.1,  3.2,  3.3, 40.4 }
        };

        public static readonly int[] globalIndices1 = { 0, 1, 2, 3 };
        public static readonly int[] globalIndices2 = { 2, 3, 4, 5 };
        public static readonly int[] globalIndices3 = { 4, 5, 6, 7 };

        public static readonly double[,] globalMatrix =
        {
            { 20.1,  1.1,  2.1,  3.1,  0.0,  0.0,  0.0,  0.0 },
            {  1.1, 20.2,  2.2,  3.2,  0.0,  0.0,  0.0,  0.0 },
            {  2.1,  2.2, 50.4,  3.3,  2.1,  3.1,  0.0,  0.0 },
            {  3.1,  3.2,  3.3, 50.6,  2.2,  3.2,  0.0,  0.0 },
            {  0.0,  0.0,  2.1,  2.2, 70.4,  1.1,  2.1,  3.1 },
            {  0.0,  0.0,  3.1,  3.2,  1.1, 70.6,  2.2,  3.2 },
            {  0.0,  0.0,  0.0,  0.0,  2.1,  2.2, 40.3,  3.3 },
            {  0.0,  0.0,  0.0,  0.0,  3.1,  3.2,  3.3, 40.4 }

        };

        public static void BuildPattern()
        {
            var pattern = SparsityPatternSymmetricColMajor.CreateEmpty(globalOrder);
            pattern.ConnectIndices(globalIndices1, true);
            pattern.ConnectIndices(globalIndices2, true);
            pattern.ConnectIndices(globalIndices3, true);

            Console.WriteLine("Sparsity pattern of the global matrix:");
            (new SparsityPatternWriter(pattern)).WriteToConsole();
        }
    }
}
