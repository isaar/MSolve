using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="SkylineMatrix"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class SkylineMatrixTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Fact]
        private static void TestArrayCopy()
        {
            var skyline = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order,
                SparsePosDef10by10.SkylineValues, SparsePosDef10by10.SkylineDiagOffsets, true, true);
            comparer.AssertEqual(SparsePosDef10by10.Matrix, skyline.CopyToArray2D());
        }

        [Fact]
        private static void TestClear()
        {
            var zero = Matrix.CreateZero(SparsePosDef10by10.Order, SparsePosDef10by10.Order);
            var skyline = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order,
                SparsePosDef10by10.SkylineValues, SparsePosDef10by10.SkylineDiagOffsets, true, true);
            skyline.Clear();
            comparer.AssertEqual(zero, skyline);
        }

        [Fact]
        private static void TestEquality()
        {
            var full = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
            var skyline = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order, 
                SparsePosDef10by10.SkylineValues, SparsePosDef10by10.SkylineDiagOffsets, true, true);
            Assert.True(skyline.Equals(full));
        }

        //TODO: Skyline operations are not part of the MKL SparseBLAS provider yet. There are only provided by the managed provider 
        //[Theory]
        //[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        [Fact]
        private static void TestMatrixVectorMultiplication(/*LinearAlgebraProviderChoice providers*/)
        {
            //TestSettings.RunMultiproviderTest(providers, delegate () {
            //

            var A = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order,
                SparsePosDef10by10.SkylineValues, SparsePosDef10by10.SkylineDiagOffsets, true, true);
            var x = Vector.CreateFromArray(SparsePosDef10by10.Lhs);
            var bExpected = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
            IVector bComputed = A.Multiply(x, false);
            comparer.AssertEqual(bExpected, bComputed);

            //});
        }
    }
}
