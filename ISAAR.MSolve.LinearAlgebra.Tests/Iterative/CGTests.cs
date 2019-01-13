using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Iterative
{
    /// <summary>
    /// Tests for <see cref="PcgAlgorithm"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class CGTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-5);

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestDenseSystem(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
                var b = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
                var xExpected = Vector.CreateFromArray(SymmPosDef10by10.Lhs);

                var builder = new CGAlgorithm.Builder();
                builder.ResidualTolerance = 1E-7;
                builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
                var cg = builder.Build();
                var xComputed = Vector.CreateZero(A.NumRows);
                CGStatistics stats = cg.Solve(A, b, xComputed, true);
                comparer.AssertEqual(xExpected, xComputed);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestSparseSystem(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
                var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
                var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

                var builder = new CGAlgorithm.Builder();
                builder.ResidualTolerance = 1E-7;
                builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
                var cg = builder.Build();
                var xComputed = Vector.CreateZero(A.NumRows);
                CGStatistics stats = cg.Solve(A, b, xComputed, true);
                comparer.AssertEqual(xExpected, xComputed);
            });
        }
    }
}
