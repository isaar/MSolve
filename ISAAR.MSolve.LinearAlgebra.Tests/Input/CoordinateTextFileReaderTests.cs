using System.IO;
using ISAAR.MSolve.LinearAlgebra.Input;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Input
{
    /// <summary>
    /// Tests for <see cref="CoordinateTextFileReader"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class CoordinateTextFileReaderTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-10);

        [Fact]
        private static void TestRandomMatrix()
        {
            // Create the random matrix and write it to a temporary file
            DokSymmetric originalMatrix = RandomUtilities.CreateRandomMatrix(1000, 0.2);
            var coordinateWriter = new CoordinateTextFileWriter();
            coordinateWriter.NumericFormat = new ExponentialFormat { NumDecimalDigits = 10 };
            string tempFile = "temp.txt";
            coordinateWriter.WriteToFile(originalMatrix, tempFile);

            // Read the temporary file and compare it with the generated matrix
            var reader = new CoordinateTextFileReader();
            DokSymmetric readMatrix = reader.ReadFileAsDokSymmetricColMajor(tempFile);
            bool success = comparer.AreEqual(originalMatrix, readMatrix);

            // Clean up
            File.Delete(tempFile);
            Assert.True(success);
        }
    }
}
