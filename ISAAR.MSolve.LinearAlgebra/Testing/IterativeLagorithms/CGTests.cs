using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.LinearSystems;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Statistics;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Testing.IterativeLagorithms
{
    public class CGTests
    {
        private delegate (Vector solution, IterativeStatistics statistics) SolveSystem(LinearSystem sys, double tol);

        public static void Run()
        {
            //CheckSystem(GetDenseSystem(), SolveWithCG, 1e-6);
            CheckSystem(GetSparseSystem(), SolveWithCG, 1e-6);
            CheckSystem(GetDenseSystem(), SolveWithPCG, 1e-6);
            CheckSystem(GetSparseSystem(), SolveWithPCG, 1e-6);
        }

        private static void CheckSystem(LinearSystem sys, SolveSystem solver, double solveTolerance)
        {
            double checkTolerance = 1e-5;
            var comparer = new Comparer(Comparer.PrintMode.Always, checkTolerance);
            (Vector solutionComputed, IterativeStatistics statistics) = solver(sys, solveTolerance);
            Console.WriteLine(statistics.ToString());
            comparer.CheckSystemSolution(sys.MatrixArray, sys.RhsArray, sys.SolutionExpected, solutionComputed.CopyToArray());
        }

        private static (Vector solution, IterativeStatistics statistics) SolveWithCG(LinearSystem sys, double tol)
        {
            var cg = new ConjugateGradient(sys.Matrix.NumColumns, tol);
            return cg.Solve(sys.Matrix, sys.Rhs);
        }

        private static (Vector solution, IterativeStatistics statistics) SolveWithPCG(LinearSystem sys, double tol)
        {
            var pcg = new PreconditionedConjugateGradient(sys.Matrix.NumColumns, tol);
            //var preconditioner = new IdentityPreconditioner(true);
            var preconditioner = new JacobiPreconditioner(sys.Matrix.GetDiagonalAsArray());
            return pcg.Solve(sys.Matrix, sys.Rhs, preconditioner);
        }

        private static LinearSystem GetDenseSystem()
        {
            return new LinearSystem
            {
                MatrixArray = SymmPositiveDefinite.matrix,
                Matrix = Matrix.CreateFromArray(SymmPositiveDefinite.matrix),
                RhsArray = SymmPositiveDefinite.rhs,
                Rhs = Vector.CreateFromArray(SymmPositiveDefinite.rhs),
                SolutionExpected = SymmPositiveDefinite.lhs
            };
        }

        private static LinearSystem GetSparseSystem()
        {
            // TODO: Add a DOK (extension) method that takes a whole matrix and only stores the matrix[i,j] != 0 entries.
            double[,] denseMatrix = SparsePositiveDefinite.matrix;
            var dok = DOKRowMajor.CreateEmpty(SparsePositiveDefinite.order, SparsePositiveDefinite.order);
            for (int j = 0; j < SparsePositiveDefinite.order; ++j)
            {
                for (int i = 0; i < SparsePositiveDefinite.order; ++i)
                {
                    if (denseMatrix[i, j] != 0) dok[i, j] = denseMatrix[i, j];
                }
            }
            return new LinearSystem
            {
                MatrixArray = SparsePositiveDefinite.matrix,
                Matrix = dok.BuildCSRMatrix(true),
                RhsArray = SparsePositiveDefinite.rhs,
                Rhs = Vector.CreateFromArray(SparsePositiveDefinite.rhs),
                SolutionExpected = SparsePositiveDefinite.lhs
            };
        }

        private class LinearSystem
        {
            internal IMatrixView Matrix { get; set; }
            internal double[,] MatrixArray { get; set; }
            internal Vector Rhs { get; set; }
            internal double[] RhsArray { get; set; }
            internal double[] SolutionExpected { get; set; }
        }
    }
}
