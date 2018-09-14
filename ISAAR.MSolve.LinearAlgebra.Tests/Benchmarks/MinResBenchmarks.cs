using System;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms.MinRes;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Benchmarks
{
    /// <summary>
    /// Benchmarks for <see cref="MinimumResidual"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class MinResBenchmarks
    {
        private static void CheckSolution(Vector solutionExpected, Vector solutionComputed)
        {
            double error = (solutionComputed - solutionExpected).Norm2() / solutionExpected.Norm2();
            Console.WriteLine("Normalized solution error = |computed - expected| / |expected| = " + error);
        }

        public static void PreconditioningIndefinite()
        {
            Console.WriteLine("Dense pos-def system WITHOUT preconditioning:");
            (Matrix A, Vector b, Vector xExpected, IPreconditioner M) = DiagonalIndefinite.BuildIndefiniteSystem(2000);
            var minres = new MinimumResidual(A.NumRows, 1e-10, 0, true, false);

            // Without preconditioning
            (Vector xSimple, MinresStatistics statsSimple) = minres.Solve(A, b);
            Console.Write(statsSimple);
            if (xSimple != null) CheckSolution(xExpected, xSimple);
            Console.WriteLine();

            // With preconditioning
            (Vector xPrec, MinresStatistics statsPrec) = minres.Solve(A, b, M);
            Console.Write(statsPrec);
            if (xPrec != null) CheckSolution(xExpected, xPrec);
            Console.WriteLine();
        }

        public static void PreconditioningPosDefDense()
        {
            Console.WriteLine("Dense pos-def system WITHOUT preconditioning:");
            var A = Matrix.CreateFromArray(SymmPosDef10by10.matrix);
            var b = Vector.CreateFromArray(SymmPosDef10by10.rhs);
            var xExpected = Vector.CreateFromArray(SymmPosDef10by10.lhs);
            var M = new JacobiPreconditioner(A.GetDiagonalAsArray());
            var minres = new MinimumResidual(A.NumRows, 1e-10, 0, true, false);

            // Without preconditioning
            (Vector xSimple, MinresStatistics statsSimple) = minres.Solve(A, b);
            Console.Write(statsSimple);
            if (xSimple != null) CheckSolution(xExpected, xSimple);
            Console.WriteLine();

            // With preconditioning
            (Vector xPrec, MinresStatistics statsPrec) = minres.Solve(A, b, M);
            Console.Write(statsPrec);
            if (xPrec != null) CheckSolution(xExpected, xPrec);
            Console.WriteLine();
        }

        public static void PreconditioningPosDefSparse()
        {
            Console.WriteLine("Assessing correctness and efficiency of preconditioned MINRES:\n");
            var A = Matrix.CreateFromArray(SparsePosDef10by10.matrix);
            var b = Vector.CreateFromArray(SparsePosDef10by10.rhs);
            var xExpected = Vector.CreateFromArray(SparsePosDef10by10.lhs);
            var M = new JacobiPreconditioner(A.GetDiagonalAsArray());
            var minres = new MinimumResidual(A.NumRows, 1e-10, 0, true, false);

            // Without preconditioning
            Console.WriteLine("Sparse pos-def system WITHOUT preconditioning:");
            (Vector xSimple, MinresStatistics statsSimple) = minres.Solve(A, b);
            Console.Write(statsSimple);
            if (xSimple != null) CheckSolution(xExpected, xSimple);
            Console.WriteLine();

            // With preconditioning
            Console.WriteLine("Sparse pos-def system WITH preconditioning:");
            (Vector xPrec, MinresStatistics statsPrec) = minres.Solve(A, b, M);
            Console.Write(statsPrec);
            if (xPrec != null) CheckSolution(xExpected, xPrec);
            Console.WriteLine();
        }
    }
}
