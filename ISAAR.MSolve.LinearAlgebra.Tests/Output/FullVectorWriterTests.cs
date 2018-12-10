using System;
using System.IO;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Output
{
    /// <summary>
    /// Tests for <see cref="FullVectorWriter"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class FullVectorWriterTests
    {
        [Fact]
        private static void TestVector1()
        {
            var vector = Vector.CreateFromArray(TestVectors.Vector1, true);
            var writer = new FullVectorWriter(false) { ArrayFormat = Array1DFormat.Brackets };
            TestWriteOperation(vector, TestVectors.FilePath, writer);
        }

        private static void TestWriteOperation(IIndexable1D vector, string referenceFile)
            => TestWriteOperation(vector, referenceFile, new FullVectorWriter(false));

        private static void TestWriteOperation(IIndexable1D vector, string referenceFile, FullVectorWriter writer)
        {
            string tempFile = Guid.NewGuid().ToString() + ".txt";
            writer.WriteToFile(vector, tempFile);
            bool success = IOUtilities.AreFilesEquivalent(referenceFile, tempFile);
            File.Delete(tempFile);
            Assert.True(success);
        }
    }
}
