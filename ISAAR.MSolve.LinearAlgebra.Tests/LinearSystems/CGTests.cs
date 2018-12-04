using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.LinearSystems
{
    /// <summary>
    /// Tests for <see cref="PcgAlgorithm"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class CGTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-5);

        [Fact]
        private static void TestDenseSystem()
        {
            var A = Matrix.CreateFromArray(SymmPosDef10by10.matrix);
            var b = Vector.CreateFromArray(SymmPosDef10by10.rhs);
            var xExpected = Vector.CreateFromArray(SymmPosDef10by10.lhs);

            var builder = new CGAlgorithm.Builder();
            builder.ResidualTolerance = 1E-7;
            builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
            var cg = builder.Build();
            var xComputed = Vector.CreateZero(A.NumRows);
            CGStatistics stats = cg.Solve(A, b, xComputed, true);
            comparer.AssertEqual(xExpected, xComputed);
        }

        [Fact]
        private static void TestSparseSystem()
        {
            var A = Matrix.CreateFromArray(SparsePosDef10by10.matrix);
            var b = Vector.CreateFromArray(SparsePosDef10by10.rhs);
            var xExpected = Vector.CreateFromArray(SparsePosDef10by10.lhs);

            var builder = new CGAlgorithm.Builder();
            builder.ResidualTolerance = 1E-7;
            builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
            var cg = builder.Build();
            var xComputed = Vector.CreateZero(A.NumRows);
            CGStatistics stats = cg.Solve(A, b, xComputed, true);
            comparer.AssertEqual(xExpected, xComputed);
        }
    }
}
