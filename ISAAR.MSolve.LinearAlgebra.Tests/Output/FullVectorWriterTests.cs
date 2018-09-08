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
        private static void TestUpperSingular()
        {
            var vector = Vector.CreateFromArray(TestVectors.vector1, true);
            var writer = new FullVectorWriter(false) { ArrayFormat = Array1DFormat.Brackets };
            TestWriteOperation(vector, TestVectors.filePath, writer);
        }

        private static void TestWriteOperation(IIndexable1D vector, string referenceFile)
            => TestWriteOperation(vector, referenceFile, new FullVectorWriter(false));

        private static void TestWriteOperation(IIndexable1D vector, string referenceFile, FullVectorWriter writer)
        {
            string tempFile = "temp.txt";
            writer.WriteToFile(vector, tempFile);
            bool success = IOUtilities.AreFilesIdentical(referenceFile, tempFile);
            File.Delete(tempFile);
            Assert.True(success);
        }
    }
}
