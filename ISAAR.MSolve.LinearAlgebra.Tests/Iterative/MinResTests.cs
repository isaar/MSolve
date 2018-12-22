using ISAAR.MSolve.LinearAlgebra.Iterative.MinimumResidual;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Iterative
{
    /// <summary>
    /// Tests for <see cref="SparsityPatternSymmetric"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class MinResTests
    {
        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestIndefiniteSystem(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var comparer = new MatrixComparer(1E-4);
                double residualTolerance = 1e-8;
                (Matrix A, Vector b, Vector xExpected, IPreconditioner M) = DiagonalIndefinite.BuildIndefiniteSystem(20);
                var minres = new MinRes(A.NumRows, residualTolerance, 0, false, false);
                (IVector xComputed, MinresStatistics stats) = minres.Solve(A, b);
                comparer.AssertEqual(xExpected, xComputed);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestPosDefDenseSystem(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var comparer = new MatrixComparer(1E-6);
                int n = SparsePosDef10by10.Order;
                var A = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
                var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
                var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);
                var minres = new MinRes(n, 1e-10, 0, false, false);
                (IVector xComputed, MinresStatistics stats) = minres.Solve(A, b);
                comparer.AssertEqual(xExpected, xComputed);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestPosDefSparseSystem(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var comparer = new MatrixComparer(1E-6);
                int n = SparsePosDef10by10.Order;
                var A = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
                var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
                var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);
                var minres = new MinRes(n, 1e-10, 0, false, false);
                (IVector xComputed, MinresStatistics stats) = minres.Solve(A, b);
                comparer.AssertEqual(xExpected, xComputed);
            });
        }
    }
}
