using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Statistics;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

// Needs Builder pattern
namespace ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms
{
    public class PreconditionedConjugateGradient
    {
        private readonly int maxIterations;
        private readonly double tolerance;

        public PreconditionedConjugateGradient(int maxIterations, double tolerance)
        {
            this.maxIterations = maxIterations;
            this.tolerance = tolerance;
        }

        public (Vector, IterativeStatistics) Solve(IMatrixView matrix, Vector rhs, IPreconditioner preconditioner)
        {
            CheckInput(matrix, rhs);
            Vector sol = Vector.CreateZero(matrix.NumColumns); // Start from x = 0
            Vector res = rhs.Copy(); // No need to do the multiplication A*0 = 0
            return SolveInternal(matrix, preconditioner, sol, res);
        }

        public (Vector, IterativeStatistics) Solve(IMatrixView matrix, Vector rhs, IPreconditioner preconditioner, 
            Vector initialGuess)
        {
            CheckInput(matrix, rhs, initialGuess);
            Vector sol = initialGuess.Copy(); // Should I copy this?
            Vector res = rhs - matrix.MultiplyRight(sol);
            return SolveInternal(matrix, preconditioner, sol, res);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="rhs"></param>
        /// <param name="sol">The solution vector initialized appropriately. It will also be returned.</param>
        /// <returns></returns>
        private (Vector, IterativeStatistics) SolveInternal(IMatrixView matrix, IPreconditioner preconditioner, 
            Vector sol, Vector res)
        {
            Vector z = preconditioner.SolveLinearSystem(res);
            Vector dir = z.Copy(); // TODO: Do I need to copy it?
            double zrDotCurrent = z * res;
            double resNormInit = res.Norm2(); // In basic CG, I could just take the sqrt(r*r), but here I have z*r.
            double resNormRatio = 1.0;
            IterativeStatistics statistics;

            for (int i = 0; i < maxIterations; ++i)
            {
                Vector matrixTimesDir = matrix.MultiplyRight(dir);
                double step = zrDotCurrent / (dir * matrixTimesDir);
                sol.AxpyIntoThis(step, dir);
                res.AxpyIntoThis(-step, matrixTimesDir);

                resNormRatio = Math.Sqrt(res.Norm2()) / resNormInit;
                if (resNormRatio < tolerance) // resNormRatio is non negative
                {
                    statistics = new IterativeStatistics
                    {
                        AlgorithmName = "Preconditioned Conjugate Gradient",
                        HasConverged = true, IterationsRequired = i + 1, NormRatio = resNormRatio
                    };
                    return (sol, statistics);
                }

                z = preconditioner.SolveLinearSystem(res); 
                double zrDotNext = z * res; //Fletcher-Reeves formula. TODO: For variable preconditioning use Polak-Ribiere
                double beta = zrDotNext / zrDotCurrent;
                dir = z.Axpy(beta, dir);
                zrDotCurrent = zrDotNext;
            }

            statistics = new IterativeStatistics
            {
                HasConverged = false, IterationsRequired = maxIterations, NormRatio = resNormRatio
            };
            return (sol, statistics);
        }

        private static void CheckInput(IMatrixView matrix, Vector rhs, Vector initialGuess)
        {
            if (matrix.NumColumns != initialGuess.Length) throw new NonMatchingDimensionsException(
                $"The matrix is {matrix.NumRows}-by-{matrix.NumColumns}),"
                + $" while the right hand side vector is {initialGuess.Length}-by-1");
            CheckInput(matrix, rhs);
        }

        private static void CheckInput(IMatrixView matrix, Vector rhs)
        {
            if (matrix.NumRows != matrix.NumColumns) throw new NonMatchingDimensionsException(
                "The matrix must be square, symmetric and positive definite.");
            Preconditions.CheckSystemSolutionDimensions(matrix, rhs);
        }
    }
}
