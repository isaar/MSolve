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
            var minres1 = new MinimumResidual(A1.NumRows, 1e-10, true, false);
            (Vector x1, MinresStatistics stats1) = minres1.Solve(A1, b1);
            Console.WriteLine(stats1);
            if (x1 != null) comparer.CheckVectorEquality(x1Expected, x1);
            Console.WriteLine();

            // Example 2
            var A2 = Matrix.CreateFromArray(SymmPositiveDefinite.matrix);
            var b2 = Vector.CreateFromArray(SymmPositiveDefinite.rhs);
            var x2Expected = Vector.CreateFromArray(SymmPositiveDefinite.lhs);
            var minres2 = new MinimumResidual(A2.NumRows, 1e-10, true, false);
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
            var minres1 = new MinimumResidual(A.NumRows, residualTolerance, true, false);
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
            var comparer = new Comparer(Comparer.PrintMode.Always, 1e-6);
            Console.WriteLine("Assessing correctness and efficiency of preconditioned MINRES:");

            // Example 1
            var A1 = Matrix.CreateFromArray(SparsePositiveDefinite.matrix);
            var b1 = Vector.CreateFromArray(SparsePositiveDefinite.rhs);
            var x1Expected = Vector.CreateFromArray(SparsePositiveDefinite.lhs);

            // Without preconditioning
            Console.WriteLine("Sparse pos-def system WITHOUT preconditioning:");
            var minres1 = new MinimumResidual(A1.NumRows, 1e-10, true, false);
            (Vector x1, MinresStatistics stats1) = minres1.Solve(A1, b1);
            Console.WriteLine(stats1);
            if (x1 != null) comparer.CheckVectorEquality(x1Expected, x1);

            // With preconditioning
            Console.WriteLine("Sparse pos-def system WITH preconditioning:");
            var minres1Prec = new PreconditionedMinimumResidual(A1.NumRows, 1e-10, true, false);
            var M1 = new JacobiPreconditioner(A1.GetDiagonalAsArray());
            (Vector x1Prec, MinresStatistics stats1Prec) = minres1Prec.Solve(A1, b1, M1);
            Console.WriteLine(stats1Prec);
            if (x1Prec != null) comparer.CheckVectorEquality(x1Expected, x1Prec);
            Console.WriteLine();



            // Example 2
            Console.WriteLine("Dense pos-def system WITHOUT preconditioning:");
            var A2 = Matrix.CreateFromArray(SymmPositiveDefinite.matrix);
            var b2 = Vector.CreateFromArray(SymmPositiveDefinite.rhs);
            var x2Expected = Vector.CreateFromArray(SymmPositiveDefinite.lhs);

            // Without preconditioning
            var minres2 = new MinimumResidual(A2.NumRows, 1e-10, true, false);
            (Vector x2, MinresStatistics stats2) = minres2.Solve(A2, b2);
            Console.WriteLine(stats2);
            if (x2 != null) comparer.CheckVectorEquality(x2Expected, x2);

            // With preconditioning
            Console.WriteLine("DENSE pos-def system WITH preconditioning:");
            var minres2Prec = new PreconditionedMinimumResidual(A2.NumRows, 1e-10, true, false);
            var M2 = new JacobiPreconditioner(A2.GetDiagonalAsArray());
            (Vector x2Prec, MinresStatistics stats2Prec) = minres2Prec.Solve(A2, b2, M2);
            Console.WriteLine(stats2Prec);
            if (x2Prec != null) comparer.CheckVectorEquality(x2Expected, x2Prec);
            Console.WriteLine();


            // Example 3
            Console.WriteLine("Diagonal indefinite system WITHOUT preconditioning:");
            (Matrix A3, Vector b3, Vector x3Expected, IPreconditioner M3) = BuildIndefiniteSystem();

            // Without preconditioning
            var minres3 = new MinimumResidual(A3.NumRows, 1e-6, true, false);
            (Vector x3, MinresStatistics stats3) = minres3.Solve(A3, b3);
            Console.WriteLine(stats3);
            if (x3 != null) comparer.CheckVectorEquality(x3Expected, x3);

            // With preconditioning
            Console.WriteLine("Diagonal indefinite system WITH preconditioning:");
            var minres3Prec = new PreconditionedMinimumResidual(A3.NumRows, 1e-6, true, false);
            (Vector x3Prec, MinresStatistics stats3Prec) = minres3Prec.Solve(A3, b3, M3);
            Console.WriteLine(stats3Prec);
            if (x3Prec != null) comparer.CheckVectorEquality(x3Expected, x3Prec);
            Console.WriteLine();
        }

        private static (Matrix A, Vector b, Vector x, IPreconditioner M) BuildIndefiniteSystem()
        {
            // Example taken from Matlab docs. Creates a diagonal matrix with half its entries being negative.
            var A = Matrix.CreateZero(40, 40);
            for (int i = 0; i < 20; ++i) A[i, i] = 20 - i;
            for (int i = 0; i < 20; ++i) A[20 + i, 20 + i] = -1 - i;

            // x = [1, 1, ... 1]
            var x = Vector.CreateWithValue(40, 1.0);

            // b[i] = A[i, i] 
            var b = Vector.CreateZero(40);
            for (int i = 0; i < 40; ++i) b[i] = A[i, i];

            // The usual Jacobi preconditioner worsens convergence. Modifications are needed.
            var positiveDiagonal = new double[40];
            for (int i = 0; i < 40; ++i) positiveDiagonal[i] = Math.Abs(A[i, i]);
            var M = new JacobiPreconditioner(positiveDiagonal);
            return (A, b, x, M);
        }
    }
}
