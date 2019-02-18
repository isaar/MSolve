using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Orthogonalization;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Orthogonalization
{
    /// <summary>
    /// Tests for <see cref="LQFactorization"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class LQFactorizationTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestFactorsLQ(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix).Transpose();
                Matrix expectedL = Matrix.CreateFromArray(RectangularFullRank10by5.LQFactorL);
                Matrix expectedQ = Matrix.CreateFromArray(RectangularFullRank10by5.LQFactorQ);

                LQFactorization factorization = A.FactorLQ();
                Matrix computedL = factorization.GetFactorL();
                Matrix computedQ = factorization.GetFactorQ();

                comparer.AssertEqual(expectedL, computedL);
                comparer.AssertEqual(expectedQ, computedQ);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestMinNormSolution(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix).Transpose();
                var b = Vector.CreateFromArray(RectangularFullRank10by5.RhsMinNorm);
                var xExpected = Vector.CreateFromArray(RectangularFullRank10by5.LhsMinNorm);

                LQFactorization factorization = A.FactorLQ();
                Vector xComputed = factorization.SolveMinNorm(b);
                comparer.AssertEqual(xExpected, xComputed);
            });
        }
    }
}
