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
    /// Tests for <see cref="LUFactorization"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class LUFactorizationTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestDeterminantInvertiblePositive(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // invertible (rank = 10) with positive det
                var A = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);
                LUFactorization factorization = A.FactorLU();
                double det = factorization.CalcDeterminant();
                comparer.AssertEqual(SquareInvertible10by10.Determinant, det);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestDeterminantInvertibleNegative(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // Switch 2 rows to make the det negative
                var A = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);
                Vector row0 = A.GetRow(0);
                Vector row9 = A.GetRow(9);
                A.SetSubrow(0, row9);
                A.SetSubrow(9, row0);
                LUFactorization factorization = A.FactorLU();
                double det = factorization.CalcDeterminant();
                comparer.AssertEqual(-SquareInvertible10by10.Determinant, det);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestDeterminantSingular1(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // singular (rank = 8)
                var A = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
                LUFactorization factorization = A.FactorLU();
                double det = factorization.CalcDeterminant();
                comparer.AssertEqual(SquareSingular10by10.Determinant, det);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestDeterminantSingular2(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // singular (rank = 9)
                var A = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.Matrix);
                LUFactorization factorization = A.FactorLU();
                double det = factorization.CalcDeterminant();
                comparer.AssertEqual(SquareSingularSingleDeficiency10by10.Determinant, det);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestFactorsLU(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // invertible (rank = 10)
                var A1 = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);
                var expectedL1 = Matrix.CreateFromArray(SquareInvertible10by10.FactorL);
                var expectedU1 = Matrix.CreateFromArray(SquareInvertible10by10.FactorU);
                LUFactorization factorization1 = A1.FactorLU();
                Matrix computedL1 = factorization1.GetFactorL();
                Matrix computedU1 = factorization1.GetFactorU();
                comparer.AssertEqual(expectedL1, computedL1);
                comparer.AssertEqual(expectedU1, computedU1);

                // singular (rank = 8)
                var A2 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
                var expectedL2 = Matrix.CreateFromArray(SquareSingular10by10.FactorL);
                var expectedU2 = Matrix.CreateFromArray(SquareSingular10by10.FactorU);
                LUFactorization factorization2 = A2.FactorLU();
                Matrix computedL2 = factorization2.GetFactorL();
                Matrix computedU2 = factorization2.GetFactorU();
                comparer.AssertEqual(expectedL2, computedL2);
                comparer.AssertEqual(expectedU2, computedU2);

                // singular (rank = 9)
                var A3 = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.Matrix);
                var expectedL3 = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.FactorL);
                var expectedU3 = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.FactorU);
                LUFactorization factorization3 = A3.FactorLU();
                Matrix computedL3 = factorization3.GetFactorL();
                Matrix computedU3 = factorization3.GetFactorU();
                comparer.AssertEqual(expectedL3, computedL3);
                comparer.AssertEqual(expectedU3, computedU3);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestInversion(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // invertible (rank = 10)
                var A1 = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);
                var inverseA1Expected = Matrix.CreateFromArray(SquareInvertible10by10.Inverse);
                LUFactorization factorization1 = A1.FactorLU();
                Matrix inverseA1Computed = factorization1.Invert(true);
                comparer.AssertEqual(inverseA1Expected, inverseA1Computed);

                // singular (rank = 8)
                var A2 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
                LUFactorization factorization2 = A2.FactorLU();
                Assert.Throws<SingularMatrixException>(() => factorization2.Invert(true));

                // singular (rank = 9)
                var A3 = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.Matrix);
                LUFactorization factorization3 = A3.FactorLU();
                Assert.Throws<SingularMatrixException>(() => factorization3.Invert(true));
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestSystemSolution(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // invertible (rank = 10)
                var A1 = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);
                var b1 = Vector.CreateFromArray(SquareInvertible10by10.Rhs);
                var x1Expected = Vector.CreateFromArray(SquareInvertible10by10.Lhs);
                LUFactorization factorization1 = A1.FactorLU();
                Vector x1Computed = factorization1.SolveLinearSystem(b1);
                comparer.AssertEqual(x1Expected, x1Computed);

                // singular (rank = 8)
                var A2 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
                var b2 = Vector.CreateFromArray(SquareSingular10by10.Rhs);
                LUFactorization factorization2 = A2.FactorLU();
                Assert.Throws<SingularMatrixException>(() => factorization2.SolveLinearSystem(b2));

                // singular (rank = 9)
                var A3 = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.Matrix);
                var b3 = Vector.CreateFromArray(SquareSingularSingleDeficiency10by10.Rhs);
                LUFactorization factorization3 = A3.FactorLU();
                Assert.Throws<SingularMatrixException>(() => factorization3.SolveLinearSystem(b3));
            });
        }
    }
}
