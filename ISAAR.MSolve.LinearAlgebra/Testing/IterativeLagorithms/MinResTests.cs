using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Statistics;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Matlab's MINRES converges faster to a better solution. Try reorthogonalized MINRES
namespace ISAAR.MSolve.LinearAlgebra.Testing.IterativeLagorithms
{
    class MinresTests
    {
        public static void CheckSolutionDefinite()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always, 1e-6);
            Console.WriteLine("Checking correctness of MINRES for pos-def system:");

            // Example 1
            var A1 = Matrix.CreateFromArray(SparsePositiveDefinite.matrix);
            var b1 = Vector.CreateFromArray(SparsePositiveDefinite.rhs);
            var x1Expected = Vector.CreateFromArray(SparsePositiveDefinite.lhs);
            var minres1 = new MinimumResidual(A1.NumRows, 1e-10, 0, true, false);
            (Vector x1, MinresStatistics stats1) = minres1.Solve(A1, b1);
            Console.WriteLine(stats1);
            if (x1 != null) comparer.CheckVectorEquality(x1Expected, x1);
            Console.WriteLine();

            // Example 2
            var A2 = Matrix.CreateFromArray(SymmPositiveDefinite.matrix);
            var b2 = Vector.CreateFromArray(SymmPositiveDefinite.rhs);
            var x2Expected = Vector.CreateFromArray(SymmPositiveDefinite.lhs);
            var minres2 = new MinimumResidual(A2.NumRows, 1e-10, 0, true, false);
            (Vector x2, MinresStatistics stats2) = minres2.Solve(A2, b2);
            Console.WriteLine(stats2);
            if (x2 != null) comparer.CheckVectorEquality(x2Expected, x2);
            Console.WriteLine();
        }

        public static void CheckSolutionIndefinite()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always, 1e-4);
            Console.WriteLine("Checking correctness of MINRES (and CG) for symmetric indefinite system:");
            double residualTolerance = 1e-6;

            (Matrix A, Vector b, Vector xExpected, IPreconditioner M) = BuildIndefiniteSystem();
            var minres1 = new MinimumResidual(A.NumRows, residualTolerance, 0, true, false);
            (Vector xMinres, MinresStatistics statsMinres) = minres1.Solve(A, b);
            Console.WriteLine(statsMinres);
            if (xMinres != null) comparer.CheckVectorEquality(xExpected, xMinres);

            var cg = new ConjugateGradient(A.NumRows, residualTolerance);
            (Vector xCG, IterativeStatistics statsCG) = cg.Solve(A, b);
            Console.WriteLine(statsCG);
            comparer.CheckVectorEquality(xExpected, xCG);
        }


        public static void AssessPreconditioning()
        {
            Console.WriteLine("Assessing correctness and efficiency of preconditioned MINRES:\n");

            // Example 1
            var A1 = Matrix.CreateFromArray(SparsePositiveDefinite.matrix);
            var b1 = Vector.CreateFromArray(SparsePositiveDefinite.rhs);
            var x1Expected = Vector.CreateFromArray(SparsePositiveDefinite.lhs);

            // Without preconditioning
            Console.WriteLine("Sparse pos-def system WITHOUT preconditioning:");
            var minres1 = new MinimumResidual(A1.NumRows, 1e-10, 0, true, false);
            (Vector x1, MinresStatistics stats1) = minres1.Solve(A1, b1);
            Console.Write(stats1);
            if (x1 != null) CheckSolution(x1Expected, x1);
            Console.WriteLine();

            // With preconditioning
            Console.WriteLine("Sparse pos-def system WITH preconditioning:");
            var minres1Prec = new MinimumResidual(A1.NumRows, 1e-10, 0, true, false);
            var M1 = new JacobiPreconditioner(A1.GetDiagonalAsArray());
            (Vector x1Prec, MinresStatistics stats1Prec) = minres1Prec.Solve(A1, b1, M1);
            Console.Write(stats1Prec);
            if (x1Prec != null) CheckSolution(x1Expected, x1Prec);
            Console.WriteLine();



            // Example 2
            Console.WriteLine("Dense pos-def system WITHOUT preconditioning:");
            var A2 = Matrix.CreateFromArray(SymmPositiveDefinite.matrix);
            var b2 = Vector.CreateFromArray(SymmPositiveDefinite.rhs);
            var x2Expected = Vector.CreateFromArray(SymmPositiveDefinite.lhs);

            // Without preconditioning
            var minres2 = new MinimumResidual(A2.NumRows, 1e-10, 0, true, false);
            (Vector x2, MinresStatistics stats2) = minres2.Solve(A2, b2);
            Console.Write(stats2);
            if (x2 != null) CheckSolution(x2Expected, x2);
            Console.WriteLine();

            // With preconditioning
            Console.WriteLine("DENSE pos-def system WITH preconditioning:");
            var minres2Prec = new MinimumResidual(A2.NumRows, 1e-10, 0, true, false);
            var M2 = new JacobiPreconditioner(A2.GetDiagonalAsArray());
            (Vector x2Prec, MinresStatistics stats2Prec) = minres2Prec.Solve(A2, b2, M2);
            Console.Write(stats2Prec);
            if (x2Prec != null) CheckSolution(x2Expected, x2Prec);
            Console.WriteLine();


            // Example 3
            Console.WriteLine("Diagonal indefinite system WITHOUT preconditioning:");
            (Matrix A3, Vector b3, Vector x3Expected, IPreconditioner M3) = BuildIndefiniteSystem();

            // Without preconditioning
            var minres3 = new MinimumResidual(A3.NumRows, 1e-6, 0, true, false); // reorthogonalization is useless for diagonal matrices
            (Vector x3, MinresStatistics stats3) = minres3.Solve(A3, b3);
            Console.Write(stats3);
            if (x3 != null) CheckSolution(x3Expected, x3);
            Console.WriteLine();

            // With preconditioning
            Console.WriteLine("Diagonal indefinite system WITH preconditioning:");
            var minres3Prec = new MinimumResidual(A3.NumRows, 1e-6, 0, true, false);
            (Vector x3Prec, MinresStatistics stats3Prec) = minres3Prec.Solve(A3, b3, M3);
            Console.Write(stats3Prec);
            if (x3Prec != null) CheckSolution(x3Expected, x3Prec);
            Console.WriteLine();
        }

        private static (Matrix A, Vector b, Vector x, IPreconditioner M) BuildIndefiniteSystem()
        {
            // Example taken from Matlab docs. Creates a diagonal matrix with half its entries being negative.
            int n = 2000; // WARNNING: must be even
            var A = Matrix.CreateZero(n, n);
            for (int i = 0; i < n/2; ++i) A[i, i] = n/2 - i;
            for (int i = 0; i < n/2; ++i) A[n/2 + i, n/2 + i] = -1 - i;

            // x = [1, 1, ... 1]
            var x = Vector.CreateWithValue(n, 1.0);

            // b[i] = A[i, i] 
            var b = Vector.CreateZero(n);
            for (int i = 0; i < n; ++i) b[i] = A[i, i];

            // The usual Jacobi preconditioner worsens convergence. Modifications are needed.
            var positiveDiagonal = new double[n];
            for (int i = 0; i < n; ++i) positiveDiagonal[i] = Math.Abs(A[i, i]);
            var M = new JacobiPreconditioner(positiveDiagonal);
            return (A, b, x, M);
        }

        private static void CheckSolution(Vector solutionExpected, Vector solutionComputed)
        {
            double error = (solutionComputed - solutionExpected).Norm2() / solutionExpected.Norm2();
            Console.WriteLine("Normalized solution error = |computed - expected| / |expected| = " + error);
        }
    }
}
