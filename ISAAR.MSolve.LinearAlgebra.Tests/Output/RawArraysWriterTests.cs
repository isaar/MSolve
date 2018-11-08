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
            var matrix = CscMatrix.CreateFromArrays(SparseRectangular10by5.numRows, SparseRectangular10by5.numCols,
                SparseRectangular10by5.cscValues, SparseRectangular10by5.cscRowIndices, SparseRectangular10by5.cscColOffsets,
                true);
            TestWriteOperation(matrix, SparseRectangular10by5.cscArraysPath);
        }

        [Fact]
        private static void TestRectangularCSR()
        {
            var matrix = CsrMatrix.CreateFromArrays(SparseRectangular10by5.numRows, SparseRectangular10by5.numCols,
                SparseRectangular10by5.csrValues, SparseRectangular10by5.csrColIndices, SparseRectangular10by5.csrRowOffsets,
                true);
            TestWriteOperation(matrix, SparseRectangular10by5.csrArraysPath);
        }

        [Fact]
        private static void TestSkylinePosDef()
        {
            var matrix = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.order, SparsePosDef10by10.skylineValues, 
                SparsePosDef10by10.skylineDiagOffsets, true, true);
            TestWriteOperation(matrix, SparsePosDef10by10.skylineArraysPath);
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
