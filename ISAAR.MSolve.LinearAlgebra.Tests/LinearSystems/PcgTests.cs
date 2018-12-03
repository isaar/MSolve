﻿using ISAAR.MSolve.LinearAlgebra.Iterative;
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
    /// Tests for <see cref="CG"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class PcgTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-5);

        [Fact]
        private static void TestPosDefDenseSystem()
        {
            int n = SymmPosDef10by10.order;
            var A = Matrix.CreateFromArray(SymmPosDef10by10.matrix);
            var b = Vector.CreateFromArray(SymmPosDef10by10.rhs);
            var xExpected = Vector.CreateFromArray(SymmPosDef10by10.lhs);

            double tol = 1E-7;
            var pcg = new PCG(new FixedMaxIterationsProvider(n), tol);
            var M = new JacobiPreconditioner(A.GetDiagonalAsArray());
            Vector xComputed = Vector.CreateZero(A.NumRows);
            CGStatistics stats = pcg.Solve(A, M, b, xComputed, true, () => Vector.CreateZero(b.Length));
            comparer.AssertEqual(xExpected, xComputed);
        }

        [Fact]
        private static void TestPosDefSparseSystem()
        {
            int n = SparsePosDef10by10.order;
            var A = Matrix.CreateFromArray(SparsePosDef10by10.matrix);
            var b = Vector.CreateFromArray(SparsePosDef10by10.rhs);
            var xExpected = Vector.CreateFromArray(SparsePosDef10by10.lhs);

            double tol = 1E-7;
            var pcg = new PCG(new FixedMaxIterationsProvider(n), tol);
            var M = new JacobiPreconditioner(A.GetDiagonalAsArray());
            Vector xComputed = Vector.CreateZero(A.NumRows);
            CGStatistics stats = pcg.Solve(A, M, b, xComputed, true, () => Vector.CreateZero(b.Length));
            comparer.AssertEqual(xExpected, xComputed);
        }

        [Fact]
        private static void TestIndefiniteSystem()
        {
            double residualTolerance = 1e-6;
            (Matrix A, Vector b, Vector xExpected, IPreconditioner M) = DiagonalIndefinite.BuildIndefiniteSystem(20);
            var cg = new CG(new PercentageMaxIterationsProvider(1.0), residualTolerance);
            Vector xComputed = Vector.CreateZero(A.NumRows);
            CGStatistics stats = cg.Solve(A, b, xComputed, true);
            Assert.False(comparer.AreEqual(xExpected, xComputed));
        }
    }
}
