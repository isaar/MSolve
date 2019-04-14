using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;

//TODO: either use the matrices in TestData or move the example matrices there.
namespace ISAAR.MSolve.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for extension methods in <see cref="MatrixNormExtensions"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class MatrixNormTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Fact]
        private static void TestNormInf()
        {
            var A = Matrix.CreateFromArray(new double[,]
            {
                { -3, 5, 7 },
                {  2, 6, 4 },
                {  0, 2, 8 }
            });
            double normExpected = 15.0;
            double normComputed = A.NormInf();
            comparer.AssertEqual(normExpected, normComputed);
        }
    }
}
