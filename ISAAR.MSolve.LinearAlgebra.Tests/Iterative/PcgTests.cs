using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
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
    public static class PcgTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-5);

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestPosDefDenseSystem(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
                var b = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
                var xExpected = Vector.CreateFromArray(SymmPosDef10by10.Lhs);

                var builder = new PcgAlgorithm.Builder();
                builder.ResidualTolerance = 1E-7;
                builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
                var pcg = builder.Build();
                var M = new JacobiPreconditioner(A.GetDiagonalAsArray());
                Vector xComputed = Vector.CreateZero(A.NumRows);
                IterativeStatistics stats = pcg.Solve(A, M, b, xComputed, true, () => Vector.CreateZero(b.Length));
                comparer.AssertEqual(xExpected, xComputed);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestPosDefSparseSystem(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
                var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
                var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

                var builder = new PcgAlgorithm.Builder();
                builder.ResidualTolerance = 1E-7;
                builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
                var pcg = builder.Build();
                var M = new JacobiPreconditioner(A.GetDiagonalAsArray());
                Vector xComputed = Vector.CreateZero(A.NumRows);
                IterativeStatistics stats = pcg.Solve(A, M, b, xComputed, true, () => Vector.CreateZero(b.Length));
                comparer.AssertEqual(xExpected, xComputed);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestIndefiniteSystem(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                (Matrix A, Vector b, Vector xExpected, IPreconditioner M) = DiagonalIndefinite.BuildIndefiniteSystem(20);
                var builder = new CGAlgorithm.Builder();
                builder.ResidualTolerance = 1E-6;
                builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
                var cg = builder.Build();
                Vector xComputed = Vector.CreateZero(A.NumRows);
                IterativeStatistics stats = cg.Solve(A, b, xComputed, true);
                Assert.False(comparer.AreEqual(xExpected, xComputed));
            });
        }
    }
}
