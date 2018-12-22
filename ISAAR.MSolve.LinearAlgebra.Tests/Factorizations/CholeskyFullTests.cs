using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Factorizations
{
    /// <summary>
    /// Tests for <see cref="CholeskyFull"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class CholeskyFullTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestDeterminant(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // positive definite
                var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
                CholeskyFull factorization = A.FactorCholesky();
                double detComputed = factorization.CalcDeterminant();
                comparer.AssertEqual(SymmPosDef10by10.Determinant, detComputed);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestInversion(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // positive definite
                var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
                var inverseAExpected = Matrix.CreateFromArray(SymmPosDef10by10.Inverse);
                CholeskyFull factorization = A.FactorCholesky();
                Matrix inverseAComputed = factorization.Invert(true);
                comparer.AssertEqual(inverseAExpected, inverseAComputed);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestFactorization(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // positive definite
                var A1 = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
                var expectedU1 = Matrix.CreateFromArray(SymmPosDef10by10.FactorU);
                CholeskyFull factorization1 = A1.FactorCholesky();
                Matrix computedU1 = factorization1.GetFactorU();
                comparer.AssertEqual(expectedU1, computedU1);

                // singular
                var A2 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
                Assert.Throws<IndefiniteMatrixException>(() => A2.FactorCholesky());
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestSystemSolution(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // positive definite
                var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
                var b = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
                var xExpected = Vector.CreateFromArray(SymmPosDef10by10.Lhs);
                CholeskyFull factorization = A.FactorCholesky();
                Vector xComputed = factorization.SolveLinearSystem(b);
                comparer.AssertEqual(xExpected, xComputed);
            });
        }
    }
}
