using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="TriangularUpper"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class TriangularUpperTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Fact]
        private static void TestArrayCopy()
        {
            // invertible
            var A1 = TriangularUpper.CreateFromArray(UpperInvertible10by10.matrix);
            comparer.AssertEqual(UpperInvertible10by10.matrix, A1.CopyToArray2D());

            // singular
            var A2 = TriangularUpper.CreateFromArray(UpperSingular10by10.matrix);
            comparer.AssertEqual(UpperSingular10by10.matrix, A2.CopyToArray2D());
        }

        [Fact]
        private static void TestMatrixVectorMultiplication()
        {
            // invertible
            var A1 = TriangularUpper.CreateFromArray(UpperInvertible10by10.matrix);
            var x1 = Vector.CreateFromArray(UpperInvertible10by10.lhs);
            var b1Expected = Vector.CreateFromArray(UpperInvertible10by10.rhs);
            Vector b1Computed = A1.MultiplyRight(x1);
            comparer.AssertEqual(b1Expected, b1Computed);

            // singular
            var A2 = TriangularUpper.CreateFromArray(UpperSingular10by10.matrix);
            var x2 = Vector.CreateFromArray(UpperSingular10by10.lhs);
            var b2Expected = Vector.CreateFromArray(UpperSingular10by10.rhs);
            Vector b2Computed = A2.MultiplyRight(x1);
            comparer.AssertEqual(b2Expected, b2Computed);
        }

        [Fact]
        private static void TestSystemSolution()
        {
            // invertible
            var A1 = TriangularUpper.CreateFromArray(UpperInvertible10by10.matrix);
            var b1 = Vector.CreateFromArray(UpperInvertible10by10.rhs);
            var x1Expected = Vector.CreateFromArray(UpperInvertible10by10.lhs);
            Vector x1Computed = A1.SolveLinearSystem(b1);
            comparer.AssertEqual(x1Expected, x1Computed);

            // singular
            var A2 = TriangularUpper.CreateFromArray(UpperSingular10by10.matrix);
            var b2 = Vector.CreateFromArray(UpperSingular10by10.rhs);
            var x2Expected = Vector.CreateFromArray(UpperSingular10by10.lhs);
            Vector x2Computed = A2.SolveLinearSystem(b2);
            Assert.False(comparer.AreEqual(x2Expected, x2Computed));

            // invertible - solve transposed (forward substitution)
            Matrix A3 = Matrix.CreateFromArray(UpperInvertible10by10.matrix).Invert().Transpose();
            Vector x3Expected = A3 * b1;
            Vector x3Computed = A1.SolveLinearSystem(b1, true);
            comparer.AssertEqual(x3Expected, x3Computed);
        }

        [Fact]
        private static void TestTransposition()
        {
            // invertible
            var A1 = TriangularUpper.CreateFromArray(UpperInvertible10by10.matrix);
            var A1TransposeExpected = MatrixOperations.Transpose(UpperInvertible10by10.matrix);
            TriangularLower A1TransposeComputed = A1.Transpose(false);
            comparer.AssertEqual(A1TransposeExpected, A1TransposeComputed.CopyToArray2D());

            // singular
            var A2 = TriangularUpper.CreateFromArray(UpperSingular10by10.matrix);
            var A2TransposeExpected = MatrixOperations.Transpose(UpperSingular10by10.matrix);
            TriangularLower A2TransposeComputed = A2.Transpose(false);
            comparer.AssertEqual(A2TransposeExpected, A2TransposeComputed.CopyToArray2D());
        }
    }
}
