using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Orthogonalization;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Orthogonalization
{
    /// <summary>
    /// Tests for <see cref="QRFactorization"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class QRFactorizationTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestEconomyFactorsQ1R1(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                int m = RectangularFullRank10by5.NumRows;
                int n = RectangularFullRank10by5.NumCols;
                var A = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
                Matrix expectedQ1 = Matrix.CreateFromArray(RectangularFullRank10by5.QRFactorQ).GetSubmatrix(0, m, 0, n);
                Matrix expectedR1 = Matrix.CreateFromArray(RectangularFullRank10by5.QRqrFactorR).GetSubmatrix(0, n, 0, n);

                QRFactorization factorization = A.FactorQR();
                Matrix computedQ1 = factorization.GetEconomyFactorQ();
                TriangularUpper computedR1 = factorization.GetEconomyFactorR();

                comparer.AssertEqual(expectedQ1, computedQ1);
                comparer.AssertEqual(expectedR1, computedR1);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestFactorsQR(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
                Matrix expectedQ = Matrix.CreateFromArray(RectangularFullRank10by5.QRFactorQ);
                Matrix expectedR = Matrix.CreateFromArray(RectangularFullRank10by5.QRqrFactorR);

                QRFactorization factorization = A.FactorQR();
                Matrix computedQ = factorization.GetFactorQ();
                Matrix computedR = factorization.GetFactorR();

                comparer.AssertEqual(expectedQ, computedQ);
                comparer.AssertEqual(expectedR, computedR);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestLeastSquaresSolution(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
                QRFactorization factorization = A.FactorQR();

                // RHS is in the column space
                var b1 = Vector.CreateFromArray(RectangularFullRank10by5.Rhs10);
                var x1Expected = Vector.CreateFromArray(RectangularFullRank10by5.Lhs5);
                Vector x1Computed = factorization.SolveLeastSquares(b1);
                comparer.AssertEqual(x1Expected, x1Computed);

                // RHS is not in the column space
                var b2 = Vector.CreateFromArray(RectangularFullRank10by5.RhsLsq);
                var x2Expected = Vector.CreateFromArray(RectangularFullRank10by5.LhsLsq);
                Vector x2Computed = factorization.SolveLeastSquares(b2);
                comparer.AssertEqual(x2Expected, x2Computed);
            });
        }
    }
}
