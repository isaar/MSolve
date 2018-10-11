using System;
using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;

// TODO: implement triangular printer for triangular, symmetric matrices
// Clean up: using writer, always delete
namespace ISAAR.MSolve.LinearAlgebra.Tests.Output
{
    /// <summary>
    /// Tests for <see cref="FullMatrixWriter"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class FullMatrixWriterTests
    {
        [Fact]
        private static void TestLowerInvertible()
        {
            var matrix = Matrix.CreateFromArray(LowerInvertible10by10.matrix);
            TestWriteOperation(matrix, LowerInvertible10by10.filePath);
        }

        [Fact]
        private static void TestLowerSingular()
        {
            var matrix = Matrix.CreateFromArray(LowerSingular10by10.matrix);
            TestWriteOperation(matrix, LowerSingular10by10.filePath);
        }

        [Fact]
        private static void TestRectangularFullRank()
        {
            var matrix = Matrix.CreateFromArray(RectangularFullRank10by5.matrix);
            TestWriteOperation(matrix, RectangularFullRank10by5.filePath);
        }

        [Fact]
        private static void TestSparsePosDef()
        {
            var matrix = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.order, SparsePosDef10by10.skylineValues,
                SparsePosDef10by10.skylineDiagOffsets, true, true);
            var writer = new FullMatrixWriter { NumericFormat = new FixedPointFormat { MaxIntegerDigits = 2 } };
            TestWriteOperation(matrix, SparsePosDef10by10.fullFormatPath, writer);
        }

        [Fact]
        private static void TestSparseRectangularCSC()
        {
            var matrix = CscMatrix.CreateFromArrays(SparseRectangular10by5.numRows, SparseRectangular10by5.numCols,
                SparseRectangular10by5.cscValues, SparseRectangular10by5.cscRowIndices, SparseRectangular10by5.cscColOffsets,
                true);
            TestWriteOperation(matrix, SparseRectangular10by5.fullFormatPath);
        }

        [Fact]
        private static void TestSparseRectangularCSR()
        {
            var matrix = CsrMatrix.CreateFromArrays(SparseRectangular10by5.numRows, SparseRectangular10by5.numCols,
                SparseRectangular10by5.csrValues, SparseRectangular10by5.csrColIndices, SparseRectangular10by5.csrRowOffsets,
                true);
            TestWriteOperation(matrix, SparseRectangular10by5.fullFormatPath);
        }

        [Fact]
        private static void TestSquareInvertible()
        {
            var matrix = Matrix.CreateFromArray(SquareInvertible10by10.matrix);
            TestWriteOperation(matrix, SquareInvertible10by10.filePath);
        }

        [Fact]
        private static void TestSquareSingular()
        {
            var matrix = Matrix.CreateFromArray(SquareSingular10by10.matrix);
            TestWriteOperation(matrix, SquareSingular10by10.filePath);
        }

        [Fact]
        private static void TestSquareSingularSingleDeficiency()
        {
            var matrix = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.matrix);
            TestWriteOperation(matrix, SquareSingularSingleDeficiency10by10.filePath);
        }

        [Fact]
        private static void TestSymmPosDef()
        {
            var matrix = Matrix.CreateFromArray(SymmPosDef10by10.matrix);
            TestWriteOperation(matrix, SymmPosDef10by10.filePath);
        }

        [Fact]
        private static void TestSymmSingular()
        {
            var matrix = Matrix.CreateFromArray(SymmSingular10by10.matrix);
            TestWriteOperation(matrix, SymmSingular10by10.filePath);
        }

        [Fact]
        private static void TestUpperInvertible()
        {
            var matrix = Matrix.CreateFromArray(UpperInvertible10by10.matrix);
            TestWriteOperation(matrix, UpperInvertible10by10.filePath);
        }

        [Fact]
        private static void TestUpperSingular()
        {
            var matrix = Matrix.CreateFromArray(UpperSingular10by10.matrix);
            TestWriteOperation(matrix, UpperSingular10by10.filePath);
        }

        private static void TestWriteOperation(IIndexable2D matrix, string referenceFile)
            => TestWriteOperation(matrix, referenceFile, new FullMatrixWriter());

        private static void TestWriteOperation(IIndexable2D matrix, string referenceFile, FullMatrixWriter writer)
        {
            string tempFile = Guid.NewGuid().ToString() + ".txt";
            writer.WriteToFile(matrix, tempFile);
            bool success = IOUtilities.AreFilesEquivalent(referenceFile, tempFile);
            File.Delete(tempFile);
            Assert.True(success);
        }
    }
}
