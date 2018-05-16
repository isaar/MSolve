using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Statistics;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms
{
    public class ConjugateGradient
    {
        private readonly int maxIterations;
        private readonly double tolerance;

        public ConjugateGradient(int maxIterations, double tolerance)
        {
            this.maxIterations = maxIterations;
            this.tolerance = tolerance;
        }

        public (Vector, IterativeStatistics) Solve(IMatrixView matrix, Vector rhs)
        {
            CheckInput(matrix, rhs);
            Vector sol = Vector.CreateZero(matrix.NumColumns); // Start from x = 0
            Vector res = rhs.Copy(); // No need to do the multiplication A*0 = 0
            return SolveInternal(matrix, sol, res);
        }

        public (Vector, IterativeStatistics) Solve(IMatrixView matrix, Vector rhs, Vector initialGuess)
        {
            CheckInput(matrix, rhs, initialGuess);
            Vector sol = initialGuess.Copy(); // Should I copy this?
            Vector res = rhs - matrix.MultiplyRight(sol);
            return SolveInternal(matrix, sol, res);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="sol">The solution vector initialized appropriately. It will also be returned.</param>
        /// <returns></returns>
        private (Vector, IterativeStatistics) SolveInternal(IMatrixView matrix, Vector sol, Vector res)
        {
            //TODO: dot = norm * norm might be faster since I need the norm anyway. Does it reduce accuracy? Needs testing;
            
            Vector dir = res.Copy();
            double resDotCurrent = res * res;
            double resNormInit = Math.Sqrt(resDotCurrent);
            double resNormRatio = 1.0;
            IterativeStatistics statistics;

            for (int i = 0; i < maxIterations; ++i)
            {
                Vector matrixTimesDir = matrix.MultiplyRight(dir);
                double step = resDotCurrent / (dir * matrixTimesDir);
                sol.AxpyIntoThis(step, dir);
                res.AxpyIntoThis(-step, matrixTimesDir);
                double resDotNext = res * res;

                resNormRatio = Math.Sqrt(resDotNext) / resNormInit;
                if (resNormRatio < tolerance) // resNormRatio is non negative
                {
                    statistics = new IterativeStatistics
                    {
                        AlgorithmName = "Conjugate Gradient",
                        HasConverged = true,
                        IterationsRequired = i + 1,
                        NormRatio = resNormRatio
                    };
                    return (sol, statistics);
                }

                double beta = resDotNext / resDotCurrent;
                dir = res.Axpy(beta, dir);
                resDotCurrent = resDotNext;
            }

            statistics = new IterativeStatistics
            {
                AlgorithmName = "Conjugate Gradient",
                HasConverged = false,
                IterationsRequired = maxIterations,
                NormRatio = resNormRatio
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
