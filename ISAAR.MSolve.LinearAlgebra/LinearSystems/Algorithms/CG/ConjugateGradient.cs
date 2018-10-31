using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: needs to throw exceptions or at least report indefinite, nonsymmetric and singular matrices.
//TODO: initialization should be done by a vector factory, instead of new Vector(..)
namespace ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms.CG
{
    /// <summary>
    /// Implements the Conjugate Gradient algorithm for solving linear systems with symmetric positive definite matrices.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ConjugateGradient
    {
        private readonly int maxIterations;
        private readonly double residualTolerance;

        /// <summary>
        /// Initializes a new instance of <see cref="ConjugateGradient"/> with the specified convergence criteria. 
        /// If any criterion is met the algorithm will temrinate.
        /// </summary>
        /// <param name="maxIterations">The maximum number of iterations before the algorithm terminates.</param>
        /// <param name="residualTolerance">If norm2(b-A*x) / norm2(b-A*x0) &lt;= <paramref name="residualTolerance"/>, the 
        ///     algorithm will terminate, where x is the current solution vector and x0 the initial guess.</param>
        public ConjugateGradient(int maxIterations, double residualTolerance)
        {
            this.maxIterations = maxIterations;
            this.residualTolerance = residualTolerance;
        }

        /// <summary>
        /// Solves the linear system A * x = b, where A = <paramref name="matrix"/> and b = <paramref name="rhsVector"/>.
        /// Initially x = 0 and then it converges to the solution.
        /// </summary>
        /// <param name="matrix">The matrix A of the linear system A * x = b. It must be symmetric positive definite.</param>
        /// <param name="rhsVector">The right hand side vector b of the linear system A * x = b. Constraints:
        ///     <paramref name="rhsVector"/>.<see cref="IIndexable1D.Length"/> == 
        ///     <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="rhsVector"/> violates the described 
        ///     constraint.</exception>
        public (IVector solution, CGStatistics stats) Solve(IMatrixView matrix, IVectorView rhsVector)
        {
            CheckInput(matrix, rhsVector);
            Vector sol = Vector.CreateZero(matrix.NumColumns); // Start from x = 0
            IVector res = rhsVector.Copy(); // No need to do the multiplication A*0 = 0
            return SolveInternal(matrix, sol, res);
        }

        /// <summary>
        /// Solves the linear system A * x = b, where A = <paramref name="matrix"/> and b = <paramref name="rhsVector"/>.
        /// Initially x = <paramref name="initialGuess"/> and then it converges to the solution.
        /// </summary>
        /// <param name="matrix">The matrix A of the linear system A * x = b. It must be symmetric positive definite.</param>
        /// <param name="rhsVector">The right hand side vector b of the linear system A * x = b. Constraints:
        ///     <paramref name="rhsVector"/>.<see cref="IIndexable1D.Length"/> == 
        ///     <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.</param>
        /// <param name="initialGuess">The vector from which to start refining the solution vector x. Constraints:
        ///     <paramref name="rhsVector"/>.<see cref="IIndexable1D.Length"/> == 
        ///     <paramref name="matrix"/>.<see cref="IIndexable2D.NumColumns"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="rhsVector"/> or 
        ///     <paramref name="initialGuess"/> violate the described constraints.</exception>
        public (IVector solution, CGStatistics stats) Solve(IMatrixView matrix, IVectorView rhsVector, IVectorView initialGuess)
        {
            CheckInput(matrix, rhsVector, initialGuess);
            IVector sol = initialGuess.Copy(); // Should I copy this?
            IVector res = rhsVector.Subtract(matrix.MultiplyRight(sol));
            return SolveInternal(matrix, sol, res);
        }

        private (IVector solution, CGStatistics stats) SolveInternal(IMatrixView matrix, IVector sol, IVector res)
        {
            //TODO: dot = norm * norm might be faster since I need the norm anyway. Does it reduce accuracy? Needs testing;
            
            IVector dir = res.Copy();
            double resDotCurrent = res.DotProduct(res);
            double resNormInit = Math.Sqrt(resDotCurrent);
            double resNormRatio = 1.0;
            CGStatistics statistics;

            for (int i = 0; i < maxIterations; ++i)
            {
                IVector matrixTimesDir = matrix.MultiplyRight(dir);
                double step = resDotCurrent / (dir.DotProduct(matrixTimesDir));
                sol.AxpyIntoThis(dir, step);
                res.AxpyIntoThis(matrixTimesDir , -step);
                double resDotNext = res.DotProduct(res);

                resNormRatio = Math.Sqrt(resDotNext) / resNormInit;
                if (resNormRatio < residualTolerance) // resNormRatio is non negative
                {
                    statistics = new CGStatistics
                    {
                        AlgorithmName = "Conjugate Gradient",
                        HasConverged = true,
                        IterationsRequired = i + 1,
                        NormRatio = resNormRatio
                    };
                    return (sol, statistics);
                }

                double beta = resDotNext / resDotCurrent;
                dir = res.Axpy(dir, beta);
                resDotCurrent = resDotNext;
            }

            statistics = new CGStatistics
            {
                AlgorithmName = "Conjugate Gradient",
                HasConverged = false,
                IterationsRequired = maxIterations,
                NormRatio = resNormRatio
            };
            return (sol, statistics);
        }

        private static void CheckInput(IMatrixView matrix, IVectorView rhs, IVectorView initialGuess)
        {
            if (matrix.NumColumns != initialGuess.Length) throw new NonMatchingDimensionsException(
                $"The matrix is {matrix.NumRows}-by-{matrix.NumColumns}),"
                + $" while the right hand side vector is {initialGuess.Length}-by-1");
            CheckInput(matrix, rhs);
        }

        private static void CheckInput(IMatrixView matrix, IVectorView rhs)
        {
            if (matrix.NumRows != matrix.NumColumns) throw new NonMatchingDimensionsException(
                "The matrix must be square, symmetric and positive definite.");
            Preconditions.CheckSystemSolutionDimensions(matrix, rhs);
        }
    }
}
