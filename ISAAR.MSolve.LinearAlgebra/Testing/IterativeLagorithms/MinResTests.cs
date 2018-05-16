using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms;
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

            (Matrix A, Vector b, Vector xExpected) = BuildIndefiniteSystem();
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

        }

        private static (Matrix A, Vector b, Vector x) BuildIndefiniteSystem()
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

            return (A, b, x);
        }
    }
}
