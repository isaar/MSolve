using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: needs to throw exceptions or at least report indefinite, nonsymmetric and singular matrices.
//TODO: initialization should be done by a vector factory, instead of new Vector(..)
namespace ISAAR.MSolve.LinearAlgebra.Iterative.Algorithms.CG
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
        /// <param name="residualTolerance">
        /// The algorithm will terminate when norm2(b-A*x) / norm2(b-A*x0) &lt;= <paramref name="residualTolerance"/>, 
        /// where x is the current solution vector and x0 the initial guess.
        /// </param>
        public ConjugateGradient(int maxIterations, double residualTolerance)
        {
            this.maxIterations = maxIterations;
            this.residualTolerance = residualTolerance;
        }

        /// <summary>
        /// Solves the linear system A * x = b, where A = <paramref name="matrix"/> and b = <paramref name="rhs"/>.
        /// Initially x = <paramref name="initialGuess"/> and then it converges to the solution.
        /// </summary>
        /// <param name="matrix">The matrix A of the linear system A * x = b. It must be symmetric positive definite.</param>
        /// <param name="rhs">
        /// The right hand side vector b of the linear system A * x = b. Constraints:
        /// <paramref name="rhs"/>.<see cref="IIndexable1D.Length"/> 
        /// == <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.
        /// </param>
        /// <param name="solution">
        /// The vector from which to start refining the solution vector x. Constraints:
        /// <paramref name="solution"/>.<see cref="IIndexable1D.Length"/>
        /// == <paramref name="matrix"/>.<see cref="IIndexable2D.NumColumns"/>.
        /// </param>
        /// <param name="initialGuessIsZero">
        /// If <paramref name="solution"/> is 0, then set <paramref name="initialGuessIsZero"/> to true to avoid performing the
        /// operation b-A*0 before starting.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="rhs"/> or <paramref name="solution"/> violate the described constraints.
        /// </exception>
        public CGStatistics Solve(IMatrixView matrix, IVectorView rhs, IVector solution, bool initialGuessIsZero) //TODO: find a better way to handle the case x0=0
        {
            Preconditions.CheckSystemSolutionDimensions(matrix, rhs);
            Preconditions.CheckMultiplicationDimensions(matrix.NumColumns, solution.Length);

            IVector res;
            if (initialGuessIsZero) res = rhs.Copy();
            else res = rhs.Subtract(matrix.MultiplyRight(solution));
            return SolveInternal(matrix, solution, res);
        }

        private CGStatistics SolveInternal(IMatrixView matrix, IVector sol, IVector res)
        {
            //TODO: dot = norm * norm might be faster since I need the norm anyway. Does it reduce accuracy? Needs testing;
            
            IVector dir = res.Copy();
            double resDotCurrent = res.DotProduct(res);
            double resNormInit = Math.Sqrt(resDotCurrent);
            double resNormRatio = 1.0;

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
                    return new CGStatistics
                    {
                        AlgorithmName = "Conjugate Gradient",
                        HasConverged = true,
                        IterationsRequired = i + 1,
                        NormRatio = resNormRatio
                    };
                }

                double beta = resDotNext / resDotCurrent;
                dir = res.Axpy(dir, beta);
                resDotCurrent = resDotNext;
            }

            return new CGStatistics
            {
                AlgorithmName = "Conjugate Gradient",
                HasConverged = false,
                IterationsRequired = maxIterations,
                NormRatio = resNormRatio
            };
        }
    }
}
