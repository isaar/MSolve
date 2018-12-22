using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="SymmetricMatrix"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class SymmetricMatrixTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Fact]
        private static void TestArrayCopy()
        {
            // positive definite
            var A1 = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
            comparer.AssertEqual(SymmPosDef10by10.Matrix, A1.CopyToArray2D());

            // singular
            var A2 = SymmetricMatrix.CreateFromArray(SymmSingular10by10.Matrix);
            comparer.AssertEqual(SymmSingular10by10.Matrix, A2.CopyToArray2D());
        }

        [Fact]
        private static void TestClear()
        {
            var zero = Matrix.CreateZero(SymmPosDef10by10.Order, SymmPosDef10by10.Order);
            var matrix = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
            matrix.Clear();
            comparer.AssertEqual(zero, matrix);
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestMatrixVectorMultiplication(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // invertible
                var A1 = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
                var x1 = Vector.CreateFromArray(SymmPosDef10by10.Lhs);
                var b1Expected = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
                Vector b1Computed = A1.Multiply(x1);
                comparer.AssertEqual(b1Expected, b1Computed);

                // singular
                var A2 = SymmetricMatrix.CreateFromArray(SymmSingular10by10.Matrix);
                var x2 = Vector.CreateFromArray(SymmSingular10by10.Lhs);
                var b2Expected = Vector.CreateFromArray(SymmSingular10by10.Rhs);
                Vector b2Computed = A2.Multiply(x1);
                comparer.AssertEqual(b2Expected, b2Computed);
            });
        }

        [Fact]
        private static void TestTransposition()
        {
            // invertible
            var A1 = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
            var A1TransposeExpected = MatrixOperations.Transpose(SymmPosDef10by10.Matrix);
            SymmetricMatrix A1TransposeComputed = A1.Transpose(false);
            comparer.AssertEqual(A1TransposeExpected, A1TransposeComputed.CopyToArray2D());

            // singular
            var A2 = SymmetricMatrix.CreateFromArray(SymmSingular10by10.Matrix);
            var A2TransposeExpected = MatrixOperations.Transpose(SymmSingular10by10.Matrix);
            SymmetricMatrix A2TransposeComputed = A2.Transpose(false);
            comparer.AssertEqual(A2TransposeExpected, A2TransposeComputed.CopyToArray2D());
        }
    }
}
