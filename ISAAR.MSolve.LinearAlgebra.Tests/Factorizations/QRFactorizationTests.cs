using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Factorizations
{
    /// <summary>
    /// Tests for <see cref="QRFactorization"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class QRFactorizationTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Fact]
        private static void TestEconomyFactorsQ1R1()
        {
            int m = RectangularFullRank10by5.numRows;
            int n = RectangularFullRank10by5.numCols;
            var A = Matrix.CreateFromArray(RectangularFullRank10by5.matrix);
            Matrix expectedQ1 = Matrix.CreateFromArray(RectangularFullRank10by5.qrFactorQ).GetSubmatrix(0, m, 0, n);
            Matrix expectedR1 = Matrix.CreateFromArray(RectangularFullRank10by5.qrFactorR).GetSubmatrix(0, n, 0, n);

            QRFactorization factorization = A.FactorQR();
            Matrix computedQ1 = factorization.GetEconomyFactorQ();
            TriangularUpper computedR1 = factorization.GetEconomyFactorR();

            comparer.AssertEqual(expectedQ1, computedQ1);
            comparer.AssertEqual(expectedR1, computedR1);
        }

        [Fact]
        private static void TestFactorsQR()
        {
            var A = Matrix.CreateFromArray(RectangularFullRank10by5.matrix);
            Matrix expectedQ = Matrix.CreateFromArray(RectangularFullRank10by5.qrFactorQ);
            Matrix expectedR = Matrix.CreateFromArray(RectangularFullRank10by5.qrFactorR);

            QRFactorization factorization = A.FactorQR();
            Matrix computedQ = factorization.GetFactorQ();
            Matrix computedR = factorization.GetFactorR();

            comparer.AssertEqual(expectedQ, computedQ);
            comparer.AssertEqual(expectedR, computedR);
        }

        [Fact]
        private static void TestLeastSquaresSolution()
        {
            var A = Matrix.CreateFromArray(RectangularFullRank10by5.matrix);
            QRFactorization factorization = A.FactorQR();

            // RHS is in the column space
            var b1 = Vector.CreateFromArray(RectangularFullRank10by5.rhs10);
            var x1Expected = Vector.CreateFromArray(RectangularFullRank10by5.lhs5);
            Vector x1Computed = factorization.SolveLeastSquares(b1);
            comparer.AssertEqual(x1Expected, x1Computed);

            // RHS is not in the column space
            var b2 = Vector.CreateFromArray(RectangularFullRank10by5.rhsLsq);
            var x2Expected = Vector.CreateFromArray(RectangularFullRank10by5.lhsLsq);
            Vector x2Computed = factorization.SolveLeastSquares(b2);
            comparer.AssertEqual(x2Expected, x2Computed);
        }
    }
}
