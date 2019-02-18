using System;
using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Output
{
    /// <summary>
    /// Tests for <see cref="RawArraysWriter"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class RawArraysWriterTests
    {
        [Fact]
        private static void TestRectangularCSC()
        {
            var matrix = CscMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CscValues, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
                true);
            TestWriteOperation(matrix, SparseRectangular10by5.CscArraysPath);
        }

        [Fact]
        private static void TestRectangularCSR()
        {
            var matrix = CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CsrValues, SparseRectangular10by5.CsrColIndices, SparseRectangular10by5.CsrRowOffsets,
                true);
            TestWriteOperation(matrix, SparseRectangular10by5.CsrArraysPath);
        }

        [Fact]
        private static void TestSkylinePosDef()
        {
            var matrix = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order, SparsePosDef10by10.SkylineValues, 
                SparsePosDef10by10.SkylineDiagOffsets, true, true);
            TestWriteOperation(matrix, SparsePosDef10by10.SkylineArraysPath);
        }

        private static void TestWriteOperation(ISparseMatrix matrix, string referenceFile)
            => TestWriteOperation(matrix, referenceFile, new RawArraysWriter(true, true));

        private static void TestWriteOperation(ISparseMatrix matrix, string referenceFile, RawArraysWriter writer)
        {
            string tempFile = Guid.NewGuid().ToString() + ".txt";
            writer.WriteToFile(matrix, tempFile);
            bool success = IOUtilities.AreFilesEquivalent(referenceFile, tempFile);
            File.Delete(tempFile);
            Assert.True(success);
        }
    }
}
