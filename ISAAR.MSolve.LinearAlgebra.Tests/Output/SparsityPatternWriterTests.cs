using System;
using System.IO;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Output
{
    /// <summary>
    /// Tests for <see cref="SparsityPatternWriter"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class SparsityPatternWriterTests
    {
        [Fact]
        private static void TestSparsePosDef()
        {
            int n = SparsePosDef10by10.Order;
            var pattern = SparsityPatternSymmetric.CreateEmpty(n);
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (SparsePosDef10by10.Matrix[i, j] != 0) pattern.AddEntry(i, j);
                }
            }
            TestWriteOperation(pattern, SparsePosDef10by10.PatternPath);
        }

        private static void TestWriteOperation(ISparsityPattern matrix, string referenceFile)
            => TestWriteOperation(matrix, referenceFile, new SparsityPatternWriter());

        private static void TestWriteOperation(ISparsityPattern pattern, string referenceFile, SparsityPatternWriter writer)
        {
            string tempFile = Guid.NewGuid().ToString() + ".txt";
            writer.WriteToFile(pattern, tempFile);
            bool success = IOUtilities.AreFilesEquivalent(referenceFile, tempFile);
            File.Delete(tempFile);
            Assert.True(success);
        }
    }
}
