using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Triangulation.SampleImplementations;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Triangulation.SampleImplementations
{
    /// <summary>
    /// Tests for the basic implementations of cholesky factorization and related operations found in 
    /// <see cref="CholeskyFactorizations"/>.
    /// AuthorsL Serafeim Bakalakos
    /// </summary>
    public static class CholeskyFactorizationTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);
        private static readonly double pivotTolerance = 1E-7;

        [Fact]
        private static void TestFactorizeFullUpper1()
        {
            // positive definite
            var A1 = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
            var expectedU1 = Matrix.CreateFromArray(SymmPosDef10by10.FactorU);

            CholeskyFactorizations.FactorizeFullUpper1(A1.NumColumns, A1, pivotTolerance);
            var computedU1 = Matrix.CreateFromArray(Conversions.FullUpperColMajorToArray2D(A1.RawData, false));
            comparer.AssertEqual(expectedU1, computedU1);
        }

        [Fact]
        private static void TestFactorizeFullUpper2()
        {
            // positive definite
            var A1 = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
            var U1Expected = Matrix.CreateFromArray(SymmPosDef10by10.FactorU);

            CholeskyFactorizations.FactorizeFullUpper2(A1.NumColumns, A1, pivotTolerance);

            var U1Computed = Matrix.CreateFromArray(Conversions.FullUpperColMajorToArray2D(A1.RawData, false));
            comparer.AssertEqual(U1Expected, U1Computed);
        }

        [Fact]
        private static void TestSystemSolutionFullUpper()
        {
            int n = SymmPosDef10by10.Order;
            var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
            var b = SymmPosDef10by10.Rhs;
            var xExpected = SymmPosDef10by10.Lhs;

            // Factorization
            CholeskyFactorizations.FactorizeFullUpper2(A.NumColumns, A, pivotTolerance);

            // Forward substitution
            var x1Computed = new double[n];
            var x2Computed = new double[n];
            CholeskyFactorizations.ForwardSubstitutionFullUpper(n, A, b, x1Computed);
            CholeskyFactorizations.ForwardSubstitutionFullUpper(n, A, b, x2Computed);

            // Back substitution - column / vector version
            CholeskyFactorizations.BackSubstitutionFullUpper1(n, A, x1Computed);
            comparer.AssertEqual(xExpected, x1Computed);

            // Back substitution - dot version
            CholeskyFactorizations.BackSubstitutionFullUpper2(n, A, x2Computed);
            comparer.AssertEqual(xExpected, x2Computed);
        }
    }
}
