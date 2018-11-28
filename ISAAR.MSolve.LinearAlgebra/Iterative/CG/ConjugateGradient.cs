using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Iterative.ResidualUpdate;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: needs to throw exceptions or at least report indefinite, nonsymmetric and singular matrices.
//TODO: initialization should be done by a vector factory, instead of new Vector(..)
namespace ISAAR.MSolve.LinearAlgebra.Iterative.CG
{
    /// <summary>
    /// Implements the Conjugate Gradient algorithm for solving linear systems with a positive definite matrix.
    /// The implementation is based on the algorithm presented in section B2 of 
    /// "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain", Jonathan Richard Shewchuk, 1994
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ConjugateGradient
    {
        private readonly MaxIterationsProvider maxIterationsProvider;
        private readonly IResidualConvergence residualConvergence = new SimpleConvergence();
        private readonly IResidualCorrection residualCorrection = new PeriodicResidualCorrection();
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
        public ConjugateGradient(MaxIterationsProvider maxIterationsProvider, double residualTolerance)
        {
            this.maxIterationsProvider = maxIterationsProvider;
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

            // r = b - A * x
            IVector residual;
            if (initialGuessIsZero) residual = rhs.Copy();
            else residual = rhs.Subtract(matrix.MultiplyRight(solution));

            return SolveInternal(matrix, rhs, solution, residual);
        }

        private CGStatistics SolveInternal(IMatrixView matrix, IVectorView rhs, IVector solution, IVector residual)
        {
            int maxIterations = maxIterationsProvider.GetMaxIterationsForMatrix(matrix);

            // d = r
            IVector direction = residual.Copy();

            // δnew = δ0 = r * r
            double dotResidualNew = residual.DotProduct(residual);

            // This is only used as output
            double normResidualInitial = Math.Sqrt(dotResidualNew);

            // Initialize the strategy objects
            residualCorrection.Initialize(matrix, rhs);
            residualConvergence.Initialize(matrix, rhs, residualTolerance, dotResidualNew);

            for (int iteration = 0; iteration < maxIterations; ++iteration)
            {
                // q = A * d
                IVector matrixTimesDirection = matrix.MultiplyRight(direction);

                // α = δnew / (d * q)
                double stepSize = dotResidualNew / (direction.DotProduct(matrixTimesDirection));

                // x = x + α * d
                solution.AxpyIntoThis(direction, stepSize);

                // δold = δnew
                double dotResidualOld = dotResidualNew;

                // Normally the residual vector is updated as: r = r - α * q. However corrections might need to be applied.
                bool isResidualCorrected = residualCorrection.UpdateResidual(iteration, solution, residual,
                    (r) => r.AxpyIntoThis(matrixTimesDirection, -stepSize));

                // δnew = r * r
                dotResidualNew = residual.DotProduct(residual);

                // At this point we can check if CG has converged and exit, thus avoiding the uneccesary operations that follow.
                // During the convergence check, it may be necessary to correct the residual vector (if it hasn't already been
                // corrected) and its dot product.
                bool hasConverged = residualConvergence.HasConverged(solution, residual, ref dotResidualNew, isResidualCorrected,
                    r => r.DotProduct(r));
                if (hasConverged)
                {
                    return new CGStatistics
                    {
                        AlgorithmName = "Conjugate Gradient",
                        HasConverged = true,
                        NumIterationsRequired = iteration + 1,
                        NormRatio = Math.Sqrt(dotResidualNew) / normResidualInitial
                    };
                }

                // β = δnew / δold
                double beta = dotResidualNew / dotResidualOld;

                // d = r + β * d
                //TODO: benchmark the two options to find out which is faster
                //direction = residual.Axpy(direction, beta); //This allocates a new vector d, copies r and GCs the existing d.
                direction.LinearCombinationIntoThis(beta, residual, 1.0); //This performs additions instead of copying and needless multiplications.
            }

            // We reached the max iterations before CG converged
            return new CGStatistics
            {
                AlgorithmName = "Conjugate Gradient",
                HasConverged = false,
                NumIterationsRequired = maxIterations,
                NormRatio = Math.Sqrt(dotResidualNew) / normResidualInitial
            };
        }
        
        private CGStatistics SolveInternalOLD(IMatrixView matrix, IVector sol, IVector res)
        {
            //TODO: dot = norm * norm might be faster since I need the norm anyway. Does it reduce accuracy? Needs testing;

            int maxIterations = maxIterationsProvider.GetMaxIterationsForMatrix(matrix);
            IVector dir = res.Copy();
            double resDotCurrent = res.DotProduct(res);
            double resNormInit = Math.Sqrt(resDotCurrent);
            double resNormRatio = 1.0;

            for (int i = 0; i < maxIterations; ++i)
            {
                IVector matrixTimesDir = matrix.MultiplyRight(dir);
                double step = resDotCurrent / (dir.DotProduct(matrixTimesDir));
                sol.AxpyIntoThis(dir, step);
                res.AxpyIntoThis(matrixTimesDir, -step);
                double resDotNext = res.DotProduct(res);

                resNormRatio = Math.Sqrt(resDotNext) / resNormInit;
                if (resNormRatio < residualTolerance) // resNormRatio is non negative
                {
                    return new CGStatistics
                    {
                        AlgorithmName = "Conjugate Gradient",
                        HasConverged = true,
                        NumIterationsRequired = i + 1,
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
                NumIterationsRequired = maxIterations,
                NormRatio = resNormRatio
            };
        }
    }
}
