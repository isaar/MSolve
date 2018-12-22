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
    /// Tests for <see cref="CholeskyPacked"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class CholeskyPackedTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestDeterminant(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // positive definite
                var A = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
                CholeskyPacked factorization = A.FactorCholesky();
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
                var A = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
                var inverseAExpected = Matrix.CreateFromArray(SymmPosDef10by10.Inverse);
                CholeskyPacked factorization = A.FactorCholesky();
                SymmetricMatrix inverseAComputed = factorization.Invert(true);
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
                var A1 = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
                var expectedU1 = Matrix.CreateFromArray(SymmPosDef10by10.FactorU);
                CholeskyPacked factorization1 = A1.FactorCholesky();
                TriangularUpper computedU1 = factorization1.GetFactorU();
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
                var A = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
                var b = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
                var xExpected = Vector.CreateFromArray(SymmPosDef10by10.Lhs);
                CholeskyPacked factorization = A.FactorCholesky();
                Vector xComputed = factorization.SolveLinearSystem(b);
                comparer.AssertEqual(xExpected, xComputed);
            });
        }
    }
}
