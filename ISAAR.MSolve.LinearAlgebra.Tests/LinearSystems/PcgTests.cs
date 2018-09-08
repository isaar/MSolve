using ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms.CG;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.LinearSystems
{
    /// <summary>
    /// Tests for <see cref="PreconditionedConjugateGradient"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class PcgTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-5);

        [Fact]
        private static void TestDenseSystem()
        {
            int n = SymmPosDef10by10.order;
            var A = Matrix.CreateFromArray(SymmPosDef10by10.matrix);
            var b = Vector.CreateFromArray(SymmPosDef10by10.rhs);
            var xExpected = Vector.CreateFromArray(SymmPosDef10by10.lhs);

            double tol = 1E-7;
            var cg = new ConjugateGradient(n, tol);
            (Vector xComputed, CGStatistics stats) = cg.Solve(A, b);
            comparer.AssertEqual(xExpected, xComputed);
        }

        [Fact]
        private static void TestSparseSystem()
        {
            int n = SparsePosDef10by10.order;
            var A = Matrix.CreateFromArray(SparsePosDef10by10.matrix);
            var b = Vector.CreateFromArray(SparsePosDef10by10.rhs);
            var xExpected = Vector.CreateFromArray(SparsePosDef10by10.lhs);

            double tol = 1E-7;
            var cg = new ConjugateGradient(n, tol);
            (Vector xComputed, CGStatistics stats) = cg.Solve(A, b);
            comparer.AssertEqual(xExpected, xComputed);
        }
    }
}
