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

        public static readonly double[] lhs5 = { 7.5025, 9.8100, 2.3352, 0.9623, 3.8458 };
        public static readonly double[] rhs5 = { 4.832870855, 46.097793477, 19.220504358, 4.022319602, 36.194372989, 5.206109486, 6.910442137, 14.814021600, 13.978989923, 27.637938654 };

        public static readonly double[] lhs10 = { 7.5025, 9.8100, 2.3352, 0.9623, 3.8458, 5.0027, 5.7026, 9.7663, 4.9286, 4.0088 };
        public static readonly double[] rhs10 = { 38.358004392, 42.386998564, 39.584157392, 117.750301553, 40.799145145 };

        public static void CheckBuilders()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var dokRowMajor = DOKRowMajor.CreateEmpty(numRows, numCols);
            var dokColMajor = DOKColMajor.CreateEmpty(numRows, numCols);
            for (int i = 0; i < numRows; ++i)
            {
                for (int j = 0; j < numCols; ++j)
                {
                    if (matrix[i, j] != 0.0)
                    {
                        dokRowMajor[i, j] = matrix[i, j];
                        dokColMajor[i, j] = matrix[i, j];
                    }
                }
            }
            Console.WriteLine();
            Console.WriteLine("Checking DOKRowMajor indexer: ");
            comparer.CheckMatrixEquality(matrix, DenseStrategies.CopyToArray2D(dokRowMajor));
            Console.WriteLine();
            Console.WriteLine("Checking DOKColMajor indexer: ");
            comparer.CheckMatrixEquality(matrix, DenseStrategies.CopyToArray2D(dokColMajor));

            Console.WriteLine();
            Console.WriteLine("Checking DOKRowMajor to CSR conversion: ");
            comparer.CheckMatrixEquality(matrix, DenseStrategies.CopyToArray2D(dokRowMajor.BuildCSRMatrix(true)));
            Console.WriteLine();
            Console.WriteLine("Checking DOKColMajor to CSC conversion: ");
            comparer.CheckMatrixEquality(matrix, DenseStrategies.CopyToArray2D(dokColMajor.BuildCSCMatrix(true)));
        }

        public static void CheckMatrixMatrixMult()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var matrix5x5 = Matrix.CreateFromArray(SquareInvertible.matrix).Slice(0, 5, 0, 5); //TODO: add a 5x5 matrix and its products
            var matrix10x10 = Matrix.CreateFromArray(SquareInvertible.matrix);
            var thisTimesMatrix5x5 = 
                MatrixOperations.MatrixTimesMatrix(matrix, matrix5x5.CopyToArray2D());
            var thisTimesTransposeMatrix5x5 = 
                MatrixOperations.MatrixTimesMatrix(matrix, matrix5x5.Transpose().CopyToArray2D());
            var transposeThisTimesMatrix10x10 = MatrixOperations.MatrixTimesMatrix(
                MatrixOperations.Transpose(matrix), matrix10x10.CopyToArray2D());
            var transposeThisTimesTransposeMatrix10x10 = MatrixOperations.MatrixTimesMatrix(
                MatrixOperations.Transpose(matrix), matrix10x10.Transpose().CopyToArray2D());
            var matrix10x10TimesThis = 
                MatrixOperations.MatrixTimesMatrix(matrix10x10.CopyToArray2D(), matrix);
            var transposeMatrix10x10TimesThis = 
                MatrixOperations.MatrixTimesMatrix(matrix10x10.Transpose().CopyToArray2D(), matrix);
            var matrix5x5TimesTransposeThis = 
                MatrixOperations.MatrixTimesMatrix(matrix5x5.CopyToArray2D(), MatrixOperations.Transpose(matrix));
            var transposeMatrix5x5TimesTransposeThis =
                MatrixOperations.MatrixTimesMatrix(matrix5x5.Transpose().CopyToArray2D(), MatrixOperations.Transpose(matrix));

            var csr = CSRMatrix.CreateFromArrays(numRows, numCols, csrValues, csrColIndices, csrRowOffsets, true);
            var csc = CSCMatrix.CreateFromArrays(numRows, numCols, cscValues, cscRowIndices, cscColOffsets, true);

            // CSR: Right multiplications
            Console.WriteLine();
            Console.WriteLine("CSR * dense: ");
            Matrix csrTimesMatrix5x5 = csr.MultiplyRight(matrix5x5, false, false);
            comparer.CheckMatrixEquality(thisTimesMatrix5x5, csrTimesMatrix5x5.CopyToArray2D());

            Console.WriteLine();
            Console.WriteLine("CSR * transpose(dense): ");
            Matrix csrTimesTransposeMatrix5x5 = csr.MultiplyRight(matrix5x5, false, true);
            comparer.CheckMatrixEquality(thisTimesTransposeMatrix5x5, csrTimesTransposeMatrix5x5.CopyToArray2D());

            Console.WriteLine();
            Console.WriteLine("transpose(CSR) * dense: ");
            Matrix transposeCsrTimesMatrix10x10 = csr.MultiplyRight(matrix10x10, true, false);
            comparer.CheckMatrixEquality(transposeThisTimesMatrix10x10, transposeCsrTimesMatrix10x10.CopyToArray2D());

            Console.WriteLine();
            Console.WriteLine("transpose(CSR) * transpose(dense): ");
            Matrix transposeCsrTimesTransposeMatrix10x10 = csr.MultiplyRight(matrix10x10, true, true);
            comparer.CheckMatrixEquality(transposeThisTimesTransposeMatrix10x10, transposeCsrTimesTransposeMatrix10x10.CopyToArray2D());

            //CSR: Left multiplications
            Console.WriteLine();
            Console.WriteLine("dense * CSR: ");
            Matrix matrix10x10TimesCSR = csr.MultiplyLeft(matrix10x10, false, false);
            comparer.CheckMatrixEquality(matrix10x10TimesThis, matrix10x10TimesCSR.CopyToArray2D());

            Console.WriteLine();
            Console.WriteLine("transpose(dense) * CSR: ");
            Matrix transposeMatrix10x10TimesCSR = csr.MultiplyLeft(matrix10x10, false, true);
            comparer.CheckMatrixEquality(transposeMatrix10x10TimesThis, transposeMatrix10x10TimesCSR.CopyToArray2D());

            Console.WriteLine();
            Console.WriteLine("dense * transpose(CSR): ");
            Matrix matrix5x5TimesTransposeCSR = csr.MultiplyLeft(matrix5x5, true, false);
            comparer.CheckMatrixEquality(matrix5x5TimesTransposeThis, matrix5x5TimesTransposeCSR.CopyToArray2D());

            Console.WriteLine();
            Console.WriteLine("transpose(dense) * transpose(CSR): ");
            Matrix transposeMatrix5x5TimesTransposeCSR = csr.MultiplyLeft(matrix5x5, true, true);
            comparer.CheckMatrixEquality(transposeMatrix5x5TimesTransposeThis, transposeMatrix5x5TimesTransposeCSR.CopyToArray2D());

            // CSC: Right multiplications
            Console.WriteLine();
            Console.WriteLine("CSC * dense: ");
            Matrix cscTimesMatrix5x5 = csc.MultiplyRight(matrix5x5, false, false);
            comparer.CheckMatrixEquality(thisTimesMatrix5x5, cscTimesMatrix5x5.CopyToArray2D());

            Console.WriteLine();
            Console.WriteLine("CSC * transpose(dense): ");
            Matrix cscTimesTransposeMatrix5x5 = csc.MultiplyRight(matrix5x5, false, true);
            comparer.CheckMatrixEquality(thisTimesTransposeMatrix5x5, cscTimesTransposeMatrix5x5.CopyToArray2D());

            Console.WriteLine();
            Console.WriteLine("transpose(CSC) * dense: ");
            Matrix transposeCscTimesMatrix10x10 = csc.MultiplyRight(matrix10x10, true, false);
            comparer.CheckMatrixEquality(transposeThisTimesMatrix10x10, transposeCscTimesMatrix10x10.CopyToArray2D());

            Console.WriteLine();
            Console.WriteLine("transpose(CSC) * transpose(dense): ");
            Matrix transposeCscTimesTransposeMatrix10x10 = csc.MultiplyRight(matrix10x10, true, true);
            comparer.CheckMatrixEquality(transposeThisTimesTransposeMatrix10x10, transposeCscTimesTransposeMatrix10x10.CopyToArray2D());

            //CSC: Left multiplications
            Console.WriteLine();
            Console.WriteLine("dense * CSC: ");
            Matrix matrix10x10TimesCSC = csc.MultiplyLeft(matrix10x10, false, false);
            comparer.CheckMatrixEquality(matrix10x10TimesThis, matrix10x10TimesCSC.CopyToArray2D());

            Console.WriteLine();
            Console.WriteLine("transpose(dense) * CSC: ");
            Matrix transposeMatrix10x10TimesCSC = csc.MultiplyLeft(matrix10x10, false, true);
            comparer.CheckMatrixEquality(transposeMatrix10x10TimesThis, transposeMatrix10x10TimesCSC.CopyToArray2D());

            Console.WriteLine();
            Console.WriteLine("dense * transpose(CSC): ");
            Matrix matrix5x5TimesTransposeCSC = csc.MultiplyLeft(matrix5x5, true, false);
            comparer.CheckMatrixEquality(matrix5x5TimesTransposeThis, matrix5x5TimesTransposeCSC.CopyToArray2D());

            Console.WriteLine();
            Console.WriteLine("transpose(dense) * transpose(CSC): ");
            Matrix transposeMatrix5x5TimesTransposeCSC = csc.MultiplyLeft(matrix5x5, true, true);
            comparer.CheckMatrixEquality(transposeMatrix5x5TimesTransposeThis, transposeMatrix5x5TimesTransposeCSC.CopyToArray2D());
        }

        public static void CheckMatrixVectorMult()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var csr = CSRMatrix.CreateFromArrays(numRows, numCols, csrValues, csrColIndices, csrRowOffsets, true);
            var csc = CSCMatrix.CreateFromArrays(numRows, numCols, cscValues, cscRowIndices, cscColOffsets, true);
            var x5 = Vector.CreateFromArray(lhs5);
            var x10 = Vector.CreateFromArray(lhs10);

            Vector csrTimesX5 = csr.MultiplyRight(x5, false);
            Vector cscTimesX5 = csc.MultiplyRight(x5, false);
            Vector csrTimesX10 = csr.MultiplyRight(x10, true);
            Vector cscTimesX10 = csc.MultiplyRight(x10, true);

            comparer.CheckMatrixVectorMult(matrix, lhs5, rhs5, csrTimesX5.InternalData);
            comparer.CheckMatrixVectorMult(matrix, lhs5, rhs5, cscTimesX5.InternalData);
            comparer.CheckMatrixVectorMult(matrix, lhs10, rhs10, csrTimesX10.InternalData);
            comparer.CheckMatrixVectorMult(matrix, lhs10, rhs10, cscTimesX10.InternalData);

            //TODO: the next should be hard-coded
            var x15 = Vector.CreateWithValue(5, 1000.0).Append(x5);
            var b20 = Vector.CreateWithValue(20, 10.0);
            var b20Expected = Vector.CreateWithValue(10, 10.0).Append(csrTimesX5);
            csr.MultiplyVectorSection(x15, 5, b20, 10);
            comparer.CheckVectorEquality(b20Expected, b20);
        }

        public static void CheckIndexing()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var csr = CSRMatrix.CreateFromArrays(numRows, numCols, csrValues, csrColIndices, csrRowOffsets, true);
            var csc = CSCMatrix.CreateFromArrays(numRows, numCols, cscValues, cscRowIndices, cscColOffsets, true);
            comparer.CheckMatrixEquality(matrix, csr.CopyToFullMatrix().CopyToArray2D());
            comparer.CheckMatrixEquality(matrix, csc.CopyToFullMatrix().CopyToArray2D());
        }

        public static void CheckEquals()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var full = Matrix.CreateFromArray(matrix);
            var csr = CSRMatrix.CreateFromArrays(numRows, numCols, csrValues, csrColIndices, csrRowOffsets, true);
            var csc = CSCMatrix.CreateFromArrays(numRows, numCols, cscValues, cscRowIndices, cscColOffsets, true);
            if (csr.Equals(full)) Console.WriteLine("CSR.Equals() works fine.");
            else Console.WriteLine("Error in CSR.Equals().");
            if (csc.Equals(full)) Console.WriteLine("CSC.Equals() works fine.");
            else Console.WriteLine("Error in CSC.Equals().");
            if ( (full.Equals(csr)) && (full.Equals(csc)) ) Console.WriteLine("Matrix.Equals() works fine.");
            else Console.WriteLine("Error in Matrix.Equals().");
        }

        public static void CheckTransposition()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var csr = CSRMatrix.CreateFromArrays(numRows, numCols, csrValues, csrColIndices, csrRowOffsets, true);
            var csc = CSCMatrix.CreateFromArrays(numRows, numCols, cscValues, cscRowIndices, cscColOffsets, true);
            var csrTrans = csr.Transpose(false);
            var cscTrans = csc.Transpose(false);
            comparer.CheckMatrixEquality(MatrixOperations.Transpose(matrix), csrTrans.CopyToFullMatrix().CopyToArray2D());
            comparer.CheckMatrixEquality(MatrixOperations.Transpose(matrix), cscTrans.CopyToFullMatrix().CopyToArray2D());
        }

        public static void Print()
        {
            var csr = CSRMatrix.CreateFromArrays(numRows, numCols, csrValues, csrColIndices, csrRowOffsets, true);
            var csc = CSCMatrix.CreateFromArrays(numRows, numCols, cscValues, cscRowIndices, cscColOffsets, true);

            Console.WriteLine("CSR (full) = ");
            (new FullMatrixWriter(csr)).WriteToConsole();
            Console.WriteLine();
            Console.WriteLine("CSR (arrays) = ");
            (new RawArraysWriter(csr, false)).WriteToConsole();
            Console.WriteLine();
            Console.WriteLine("CSR (sparse entries) = ");
            (new CoordinateTextFileWriter(csr)).WriteToConsole();
            Console.WriteLine();

            Console.WriteLine("CSC (full) = ");
            (new FullMatrixWriter(csc)).WriteToConsole();
            Console.WriteLine();
            Console.WriteLine("CSC (arrays) = ");
            (new RawArraysWriter(csc, false)).WriteToConsole();
            Console.WriteLine();
            Console.WriteLine("CSC (sparse entries) = ");
            (new CoordinateTextFileWriter(csc)).WriteToConsole();
            Console.WriteLine();
        }
    }
}
