using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices
{
    public class SparsePositiveDefinite
    {
        public const int order = 10;

        public static readonly double[,] matrix = {
            { 21.0,  1.0,  0.0,  4.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
            {  1.0, 22.0,  2.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0 },
            {  0.0,  2.0, 23.0,  1.0,  3.0,  1.0,  0.0,  1.0,  0.0,  0.0 },
            {  4.0,  0.0,  1.0, 24.0,  2.0,  4.0,  0.0,  0.0,  0.0,  0.0 },
            {  0.0,  0.0,  3.0,  2.0, 25.0,  5.0,  2.0,  0.0,  0.0,  1.0 },
            {  0.0,  0.0,  1.0,  4.0,  5.0, 26.0,  0.0,  0.0,  2.0,  3.0 },
            {  0.0,  1.0,  0.0,  0.0,  2.0,  0.0, 27.0,  3.0,  0.0,  0.0 },
            {  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  3.0, 28.0,  4.0,  2.0 },
            {  0.0,  0.0,  0.0,  0.0,  0.0,  2.0,  0.0,  4.0, 29.0,  0.0 },
            {  0.0,  0.0,  0.0,  0.0,  1.0,  3.0,  0.0,  2.0,  0.0, 30.0 }};

        public static readonly double[] skylineValues = {
            21.0, 22.0, 1.0,
            23.0, 2.0,
            24.0, 1.0, 0.0, 4.0,
            25.0, 2.0, 3.0,
            26.0, 5.0, 4.0, 1.0,
            27.0, 0.0, 2.0, 0.0, 0.0, 1.0,
            28.0, 3.0, 0.0, 0.0, 0.0, 1.0,
            29.0, 4.0, 0.0, 2.0,
            30.0, 0.0, 2.0, 0.0, 3.0, 1.0 };

        public static readonly int[] skylineDiagOffsets = { 0, 1, 3, 5, 9, 12, 16, 22, 28, 32, 38 };

        public static readonly double[] lhs = { 7.5025, 9.8100, 2.3352, 0.9623, 3.8458, 5.0027, 5.7026, 9.7663, 4.9286, 4.0088 };
        public static readonly double[] rhs = {
            171.2117, 233.6955, 100.5983, 83.1428, 145.5027, 177.3672, 200.7707, 320.6314, 192.0000, 158.6505 };

        public static readonly int[] matlabPermutationAMD = { 0, 1, 8, 9, 7, 2, 3, 4, 5, 6 };

        public static void Print()
        {
            var skyline = SkylineMatrix.CreateFromArrays(order, skylineValues, skylineDiagOffsets, true, true);
            var fullWriter = new FullMatrixWriter { NumericFormat = new FixedPointFormat { MaxIntegerDigits = 2 } };
            Console.WriteLine("Skyline (full) = ");
            fullWriter.WriteToConsole(skyline);
            Console.WriteLine();

            Console.WriteLine("Skyline (arrays) = ");
            (new RawArraysWriter(false)).WriteToConsole(skyline);
            Console.WriteLine();

            Console.WriteLine("Skyline (sparse entries) = ");
            (new CoordinateTextFileWriter()).WriteToConsole(skyline);
            Console.WriteLine();
        }

        public static void PrintPatternAsBoolean()
        {
            var pattern = SparsityPatternSymmetric.CreateEmpty(order);
            for (int i = 0; i < order; ++i)
            {
                for (int j = 0; j < order; ++j)
                {
                    if (matrix[i, j] != 0) pattern.AddEntry(i, j);
                }
            }
            Console.WriteLine("Sparsity pattern of the matrix:");
            (new SparsityPatternWriter()).WriteToConsole(pattern);
        }
    }
}
