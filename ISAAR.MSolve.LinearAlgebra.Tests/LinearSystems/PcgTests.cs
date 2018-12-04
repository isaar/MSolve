using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.LinearSystems
{
    /// <summary>
    /// Tests for <see cref="CGAlgorithm"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class PcgTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-5);

        [Fact]
        private static void TestPosDefDenseSystem()
        {
            var A = Matrix.CreateFromArray(SymmPosDef10by10.matrix);
            var b = Vector.CreateFromArray(SymmPosDef10by10.rhs);
            var xExpected = Vector.CreateFromArray(SymmPosDef10by10.lhs);

            var builder = new PcgAlgorithm.Builder();
            builder.ResidualTolerance = 1E-7;
            builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
            var pcg = builder.Build();
            var M = new JacobiPreconditioner(A.GetDiagonalAsArray());
            Vector xComputed = Vector.CreateZero(A.NumRows);
            CGStatistics stats = pcg.Solve(A, M, b, xComputed, true, () => Vector.CreateZero(b.Length));
            comparer.AssertEqual(xExpected, xComputed);
        }

        [Fact]
        private static void TestPosDefSparseSystem()
        {
            var A = Matrix.CreateFromArray(SparsePosDef10by10.matrix);
            var b = Vector.CreateFromArray(SparsePosDef10by10.rhs);
            var xExpected = Vector.CreateFromArray(SparsePosDef10by10.lhs);

            var builder = new PcgAlgorithm.Builder();
            builder.ResidualTolerance = 1E-7;
            builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
            var pcg = builder.Build();
            var M = new JacobiPreconditioner(A.GetDiagonalAsArray());
            Vector xComputed = Vector.CreateZero(A.NumRows);
            CGStatistics stats = pcg.Solve(A, M, b, xComputed, true, () => Vector.CreateZero(b.Length));
            comparer.AssertEqual(xExpected, xComputed);
        }

        [Fact]
        private static void TestIndefiniteSystem()
        {
            (Matrix A, Vector b, Vector xExpected, IPreconditioner M) = DiagonalIndefinite.BuildIndefiniteSystem(20);
            var builder = new CGAlgorithm.Builder();
            builder.ResidualTolerance = 1E-6;
            builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
            var cg = builder.Build();
            Vector xComputed = Vector.CreateZero(A.NumRows);
            CGStatistics stats = cg.Solve(A, b, xComputed, true);
            Assert.False(comparer.AreEqual(xExpected, xComputed));
        }
    }
}
