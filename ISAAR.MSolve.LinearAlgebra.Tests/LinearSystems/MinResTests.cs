using ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms.MinRes;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.LinearSystems
{
    /// <summary>
    /// Tests for <see cref="SparsityPatternSymmetric"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class MinResTests
    {
        [Fact]
        private static void TestIndefiniteSystem()
        {
            var comparer = new MatrixComparer(1E-4);
            double residualTolerance = 1e-8;
            (Matrix A, Vector b, Vector xExpected, IPreconditioner M) = DiagonalIndefinite.BuildIndefiniteSystem(20);
            var minres = new MinimumResidual(A.NumRows, residualTolerance, 0, false, false);
            (Vector xComputed, MinresStatistics stats) = minres.Solve(A, b);
            comparer.AssertEqual(xExpected, xComputed);
        }

        [Fact]
        private static void TestPosDefDenseSystem()
        {
            var comparer = new MatrixComparer(1E-6);
            int n = SparsePosDef10by10.order;
            var A = Matrix.CreateFromArray(SparsePosDef10by10.matrix);
            var b = Vector.CreateFromArray(SparsePosDef10by10.rhs);
            var xExpected = Vector.CreateFromArray(SparsePosDef10by10.lhs);
            var minres = new MinimumResidual(n, 1e-10, 0, false, false);
            (Vector xComputed, MinresStatistics stats) = minres.Solve(A, b);
            comparer.AssertEqual(xExpected, xComputed);
        }

        [Fact]
        private static void TestPosDefSparseSystem()
        {
            var comparer = new MatrixComparer(1E-6);
            int n = SparsePosDef10by10.order;
            var A = Matrix.CreateFromArray(SparsePosDef10by10.matrix);
            var b = Vector.CreateFromArray(SparsePosDef10by10.rhs);
            var xExpected = Vector.CreateFromArray(SparsePosDef10by10.lhs);
            var minres = new MinimumResidual(n, 1e-10, 0, false, false);
            (Vector xComputed, MinresStatistics stats) = minres.Solve(A, b);
            comparer.AssertEqual(xExpected, xComputed);

        }
    }
}
