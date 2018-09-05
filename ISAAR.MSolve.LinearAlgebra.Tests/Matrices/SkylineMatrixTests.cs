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
            var skyline = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.order,
                SparsePosDef10by10.skylineValues, SparsePosDef10by10.skylineDiagOffsets, true, true);
            comparer.AssertEqual(SparsePosDef10by10.matrix, skyline.CopyToArray2D());
        }

        [Fact]
        private static void TestEquality()
        {
            var full = Matrix.CreateFromArray(SparsePosDef10by10.matrix);
            var skyline = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.order, 
                SparsePosDef10by10.skylineValues, SparsePosDef10by10.skylineDiagOffsets, true, true);
            Assert.True(skyline.Equals(full));
        }

        [Fact]
        private static void TestMatrixVectorMultiplication()
        {
            var A = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.order,
                SparsePosDef10by10.skylineValues, SparsePosDef10by10.skylineDiagOffsets, true, true);
            var x = Vector.CreateFromArray(SparsePosDef10by10.lhs);
            var bExpected = Vector.CreateFromArray(SparsePosDef10by10.rhs);
            Vector bComputed = A.MultiplyRight(x, false);
            comparer.AssertEqual(bExpected, bComputed);
        }
    }
}
