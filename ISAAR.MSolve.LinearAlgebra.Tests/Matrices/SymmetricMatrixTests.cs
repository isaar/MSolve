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
            var A1 = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.matrix);
            comparer.AssertEqual(SymmPosDef10by10.matrix, A1.CopyToArray2D());

            // singular
            var A2 = SymmetricMatrix.CreateFromArray(SymmSingular10by10.matrix);
            comparer.AssertEqual(SymmSingular10by10.matrix, A2.CopyToArray2D());
        }

        [Fact]
        private static void TestMatrixVectorMultiplication()
        {
            // invertible
            var A1 = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.matrix);
            var x1 = Vector.CreateFromArray(SymmPosDef10by10.lhs);
            var b1Expected = Vector.CreateFromArray(SymmPosDef10by10.rhs);
            Vector b1Computed = A1.MultiplyRight(x1);
            comparer.AssertEqual(b1Expected, b1Computed);

            // singular
            var A2 = SymmetricMatrix.CreateFromArray(SymmSingular10by10.matrix);
            var x2 = Vector.CreateFromArray(SymmSingular10by10.lhs);
            var b2Expected = Vector.CreateFromArray(SymmSingular10by10.rhs);
            Vector b2Computed = A2.MultiplyRight(x1);
            comparer.AssertEqual(b2Expected, b2Computed);
        }

        [Fact]
        private static void TestTransposition()
        {
            // invertible
            var A1 = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.matrix);
            var A1TransposeExpected = MatrixOperations.Transpose(SymmPosDef10by10.matrix);
            SymmetricMatrix A1TransposeComputed = A1.Transpose(false);
            comparer.AssertEqual(A1TransposeExpected, A1TransposeComputed.CopyToArray2D());

            // singular
            var A2 = SymmetricMatrix.CreateFromArray(SymmSingular10by10.matrix);
            var A2TransposeExpected = MatrixOperations.Transpose(SymmSingular10by10.matrix);
            SymmetricMatrix A2TransposeComputed = A2.Transpose(false);
            comparer.AssertEqual(A2TransposeExpected, A2TransposeComputed.CopyToArray2D());
        }
    }
}
