using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Factorizations
{
    /// <summary>
    /// Tests for <see cref="LQFactorization"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class LQFactorizationTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Fact]
        private static void TestFactorsLQ()
        {
            var A = Matrix.CreateFromArray(RectangularFullRank10by5.matrix).Transpose();
            Matrix expectedL = Matrix.CreateFromArray(RectangularFullRank10by5.lqFactorL);
            Matrix expectedQ = Matrix.CreateFromArray(RectangularFullRank10by5.lqFactorQ);

            LQFactorization factorization = A.FactorLQ();
            Matrix computedL = factorization.GetFactorL();
            Matrix computedQ = factorization.GetFactorQ();

            comparer.AssertEqual(expectedL, computedL);
            comparer.AssertEqual(expectedQ, computedQ);
        }

        [Fact]
        private static void TestMinNormSolution()
        {
            var A = Matrix.CreateFromArray(RectangularFullRank10by5.matrix).Transpose();
            var b = Vector.CreateFromArray(RectangularFullRank10by5.rhsMinNorm);
            var xExpected = Vector.CreateFromArray(RectangularFullRank10by5.lhsMinNorm);

            LQFactorization factorization = A.FactorLQ();
            Vector xComputed = factorization.SolveMinNorm(b);
            comparer.AssertEqual(xExpected, xComputed);
        }
    }
}
