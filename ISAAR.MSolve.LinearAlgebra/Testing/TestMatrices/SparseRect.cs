using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices
{
    class SparseRect
    {
        public static readonly int numRows = 10;
        public static readonly int numCols = 5;

        public static readonly double[,] matrix = {
            { 0.00000,   0.00000,   1.20805,   2.09065,   0.00000 },
            { 0.72602,   3.72194,   0.64231,   2.74205,   0.00000 },
            { 1.68190,   0.00000,   2.78180,   0.00000,   0.02756 },
            { 0.00000,   0.00000,   0.83011,   0.00000,   0.54185 },
            { 2.44372,   1.52758,   0.00000,   2.98743,   0.00000 },
            { 0.00000,   0.00000,   0.77108,   3.53890,   0.00000 },
            { 0.53085,   0.00000,   0.00000,   3.04244,   0.00000 },
            { 0.00000,   0.00000,   0.00000,   0.00000,   3.85200 },
            { 0.42069,   0.00000,   2.65133,   2.70976,   0.52623 },
            { 3.19534,   0.00000,   0.00000,   3.80848,   0.00000 }};

        public static readonly double[] csrValues = {
            1.20805, 2.09065,
            0.72602, 3.72194, 0.64231, 2.74205,
            1.68190, 2.78180, 0.02756,
            0.83011, 0.54185,
            2.44372, 1.52758, 2.98743,
            0.77108, 3.53890,
            0.53085, 3.04244,
            3.85200,
            0.42069, 2.65133, 2.70976, 0.52623,
            3.19534, 3.80848 };

        public static readonly int[] csrColIndices = {
            2, 3,
            0, 1, 2, 3,
            0, 2, 4,
            2, 4,
            0, 1, 3,
            2, 3,
            0, 3,
            4,
            0, 2, 3, 4,
            0, 3 };

        public static readonly int[] csrRowOffsets = { 0, 2, 6, 9, 11, 14, 16, 18, 19, 23, 25 };


        public static readonly double[] cscValues = {
            0.72602, 1.68190, 2.44372, 0.53085, 0.42069, 3.19534,
            3.72194, 1.52758,
            1.20805, 0.64231, 2.78180, 0.83011, 0.77108, 2.65133,
            2.09065, 2.74205, 2.98743, 3.53890, 3.04244, 2.70976, 3.80848,
            0.02756, 0.54185, 3.85200, 0.52623 };

        public static readonly int[] cscRowIndices = {
            1, 2, 4, 6, 8, 9,
            1, 4,
            0, 1, 2, 3, 5, 8,
            0, 1, 4, 5, 6, 8, 9,
            2, 3, 7, 8 };

        public static readonly int[] cscColOffsets = { 0, 6, 8, 14, 21, 25 };


        public static void Print()
        {
            var csr = CsrMatrix.CreateFromArrays(numRows, numCols, csrValues, csrColIndices, csrRowOffsets, true);
            var csc = CscMatrix.CreateFromArrays(numRows, numCols, cscValues, cscRowIndices, cscColOffsets, true);

            var fullWriter = new FullMatrixWriter();
            var coordinateWriter = new CoordinateTextFileWriter();
            var rawWriter = new RawArraysWriter(false);

            Console.WriteLine("CSR (full) = ");
            fullWriter.WriteToConsole(csr);
            Console.WriteLine();
            Console.WriteLine("CSR (arrays) = ");
            rawWriter.WriteToConsole(csr);
            Console.WriteLine();
            Console.WriteLine("CSR (sparse entries) = ");
            coordinateWriter.WriteToConsole(csr);
            Console.WriteLine();

            Console.WriteLine("CSC (full) = ");
            fullWriter.WriteToConsole(csc);
            Console.WriteLine();
            Console.WriteLine("CSC (arrays) = ");
            rawWriter.WriteToConsole(csc);
            Console.WriteLine();
            Console.WriteLine("CSC (sparse entries) = ");
            coordinateWriter.WriteToConsole(csc);
            Console.WriteLine();
        }
    }
}
