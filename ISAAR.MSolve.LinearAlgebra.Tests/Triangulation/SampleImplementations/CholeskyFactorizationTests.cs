using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Triangulation.SampleImplementations;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Triangulation.SampleImplementations
{
    public static class CholeskyFactorizationTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);
        private static readonly double pivotTolerance = 1E-7;

        [Fact]
        private static void TestFullUpper1()
        {
            // positive definite
            var A1 = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
            var expectedU1 = Matrix.CreateFromArray(SymmPosDef10by10.FactorU);

            int errorCode = CholeskyFactorizations.FullUpper1(A1.NumColumns, A1, pivotTolerance);
            Assert.True(errorCode == 0);

            var computedU1 = Matrix.CreateFromArray(Conversions.FullUpperColMajorToArray2D(A1.RawData, false));
            comparer.AssertEqual(expectedU1, computedU1);
        }

        [Fact]
        private static void TestFullUpper2()
        {
            // positive definite
            var A1 = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
            var U1Expected = Matrix.CreateFromArray(SymmPosDef10by10.FactorU);

            int errorCode = CholeskyFactorizations.FullUpper2(A1.NumColumns, A1, pivotTolerance);
            Assert.True(errorCode == 0);

            var U1Computed = Matrix.CreateFromArray(Conversions.FullUpperColMajorToArray2D(A1.RawData, false));
            comparer.AssertEqual(U1Expected, U1Computed);
        }

        [Fact]
        private static void TestSystemSolutionFullUpper()
        {
            // positive definite
            int n = SymmPosDef10by10.Order;
            var A1 = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
            var b1 = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
            var x1Expected = Vector.CreateFromArray(SymmPosDef10by10.Lhs);

            int errorCode = CholeskyFactorizations.FullUpper2(A1.NumColumns, A1, pivotTolerance);
            Assert.True(errorCode == 0);

            var x1Computed = Vector.CreateZero(n);
            CholeskyFactorizations.SolveSystemFullUpper(n, A1, b1.RawData, x1Computed.RawData);

            comparer.AssertEqual(x1Expected, x1Computed);
        }
    }
}
