using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Iterative.ResidualUpdate;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Needs Builder pattern
//TODO: initialization should be done by a vector factory, instead of new Vector(..)
namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    /// <summary>
    /// Implements the Conjugate Gradient algorithm with preconditioning for solving linear systems with symmetric positive 
    /// definite matrices.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PCG
    {
        private readonly MaxIterationsProvider maxIterationsProvider;
        private readonly double residualTolerance;
        private readonly IResidualConvergence residualConvergence = new SimpleConvergence();
        private readonly IResidualCorrection residualCorrection = new NoResidualCorrection();

        /// <summary>
        /// Initializes a new instance of <see cref="PCG"/> with the specified convergence criteria. 
        /// If any criterion is met the algorithm will temrinate.
        /// </summary>
        /// <param name="maxIterations">The maximum number of iterations before the algorithm terminates.</param>
        /// <param name="residualTolerance">
        /// The algorithm will terminate when norm2(b-A*x) / norm2(b-A*x0) &lt;= <paramref name="residualTolerance"/>, 
        /// where x is the current solution vector and x0 the initial guess.
        /// </param>
        public PCG(MaxIterationsProvider maxIterationsProvider, double residualTolerance)
        {
            this.maxIterationsProvider = maxIterationsProvider;
            this.residualTolerance = residualTolerance;
        }

        /// <summary>
        /// Solves the linear system A * x = b by solving the preconditioned system inv(P) * A * inv(P)^T * y = inv(P) * b, 
        /// where A = <paramref name="matrix"/>, b = <paramref name="rhsVector"/>, x is the solution, y = P^T * x,
        /// P*P^T = <paramref name="preconditioner"/>.
        /// Initially x = <paramref name="initialGuess"/> and then it converges to the solution.
        /// </summary>
        /// <param name="matrix">The matrix A of the linear system A * x = b. It must be symmetric positive definite.</param>
        /// <param name="rhs">
        /// The right hand side vector b of the linear system A * x = b. Constraints:
        /// <paramref name="rhs"/>.<see cref="IIndexable1D.Length"/> 
        /// == <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.
        /// </param>
        /// <param name="preconditioner">
        /// A preconditioner matrix that is also symmetric positive definite and has the same dimensions as A.
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
        public CGStatistics Solve(IMatrixView matrix, IPreconditioner preconditioner, IVectorView rhs,  IVector solution,
            bool initialGuessIsZero, Func<IVector> zeroVectorInitializer) //TODO: find a better way to handle the case x0=0
        {
            Preconditions.CheckSystemSolutionDimensions(matrix, rhs);
            Preconditions.CheckMultiplicationDimensions(matrix.NumColumns, solution.Length);

            IVector res;
            if (initialGuessIsZero) res = rhs.Copy();
            else res = rhs.Subtract(matrix.MultiplyRight(solution));
            return SolveInternal(matrix, preconditioner, rhs, solution, res, zeroVectorInitializer);
        }

        private CGStatistics SolveInternal(IMatrixView matrix, IPreconditioner preconditioner, IVectorView rhs, 
            IVector solution, IVector residual, Func<IVector> zeroVectorInitializer)
        {
            int maxIterations = maxIterationsProvider.GetMaxIterationsForMatrix(matrix);

            // In contrast to the source algorithm, we initialize s here. At each iteration it will be overwritten, 
            // thus avoiding allocating deallocating a new vector.
            IVector preconditionedResidual = zeroVectorInitializer();

            // d = inv(M) * r
            IVector direction = zeroVectorInitializer();
            preconditioner.SolveLinearSystem(residual, direction);

            // δnew = δ0 = r * d
            double dotPreconditionedResidualNew = residual.DotProduct(direction);

            // This is only used as output
            double normResidualInitial = Math.Sqrt(dotPreconditionedResidualNew);

            // Initialize the strategy objects
            residualCorrection.Initialize(matrix, rhs);
            residualConvergence.Initialize(matrix, rhs, residualTolerance, dotPreconditionedResidualNew);

            for (int iteration = 0; iteration < maxIterations; ++iteration)
            {
                // q = A * d
                IVector matrixTimesDirection = matrix.MultiplyRight(direction);

                // α = δnew / (d * q)
                double stepSize = dotPreconditionedResidualNew / (direction.DotProduct(matrixTimesDirection));

                // x = x + α * d
                solution.AxpyIntoThis(direction, stepSize);

                // δold = δnew
                double dotPreconditionedResidualOld = dotPreconditionedResidualNew;

                // Normally the residual vector is updated as: r = r - α * q. However corrections might need to be applied.
                bool isResidualCorrected = residualCorrection.UpdateResidual(iteration, solution, residual,
                    (r) => r.AxpyIntoThis(matrixTimesDirection, -stepSize));

                // s = inv(M) * r
                preconditioner.SolveLinearSystem(residual, preconditionedResidual);

                // Fletcher-Reeves formula: δnew = r * s 
                //TODO: For variable preconditioning use Polak-Ribiere
                dotPreconditionedResidualNew = residual.DotProduct(preconditionedResidual);

                // At this point we can check if CG has converged and exit, thus avoiding the uneccesary operations that follow.
                // During the convergence check, it may be necessary to correct the residual vector (if it hasn't already been
                // corrected) and its dot product.
                bool hasConverged = residualConvergence.HasConverged(solution, residual, ref dotPreconditionedResidualNew,
                    isResidualCorrected, r => r.DotProduct(preconditionedResidual));
                if (hasConverged)
                {
                    return new CGStatistics
                    {
                        AlgorithmName = "Conjugate Gradient",
                        HasConverged = true,
                        NumIterationsRequired = iteration + 1,
                        NormRatio = Math.Sqrt(dotPreconditionedResidualNew) / normResidualInitial
                    };
                }

                // β = δnew / δold
                double beta = dotPreconditionedResidualNew / dotPreconditionedResidualOld;

                // d = s + β * d
                //TODO: benchmark the two options to find out which is faster
                //direction = preconditionedResidual.Axpy(direction, beta); //This allocates a new vector d, copies r and GCs the existing d.
                direction.LinearCombinationIntoThis(beta, preconditionedResidual, 1.0); //This performs additions instead of copying and needless multiplications.
            }

            // We reached the max iterations before CG converged
            return new CGStatistics
            {
                AlgorithmName = "Conjugate Gradient",
                HasConverged = false,
                NumIterationsRequired = maxIterations,
                NormRatio = Math.Sqrt(dotPreconditionedResidualNew) / normResidualInitial
            };
        }

        private CGStatistics SolveInternalOLD(IMatrixView matrix, IPreconditioner preconditioner, IVector sol, IVector res)
        {
            int maxIterations = maxIterationsProvider.GetMaxIterationsForMatrix(matrix);
            IVector z = res.Copy();
            preconditioner.SolveLinearSystem(res, z);
            IVector dir = z.Copy(); // TODO: Do I need to copy it?
            double zrDotCurrent = z.DotProduct(res);
            double resNormInit = res.Norm2(); // In basic CG, I could just take the sqrt(r*r), but here I have z*r.
            double resNormRatio = 1.0;

            for (int i = 0; i < maxIterations; ++i)
            {
                IVector matrixTimesDir = matrix.MultiplyRight(dir);
                double step = zrDotCurrent / (dir.DotProduct(matrixTimesDir));
                sol.AxpyIntoThis(dir, step);
                res.AxpyIntoThis(matrixTimesDir, - step);

                resNormRatio = Math.Sqrt(res.Norm2()) / resNormInit;
                if (resNormRatio < residualTolerance) // resNormRatio is non negative
                {
                    return new CGStatistics
                    {
                        AlgorithmName = "Preconditioned Conjugate Gradient",
                        HasConverged = true, NumIterationsRequired = i + 1, NormRatio = resNormRatio
                    };
                }

                preconditioner.SolveLinearSystem(res, z); 
                double zrDotNext = z.DotProduct(res); //Fletcher-Reeves formula. TODO: For variable preconditioning use Polak-Ribiere
                double beta = zrDotNext / zrDotCurrent;
                dir = z.Axpy(dir, beta);
                zrDotCurrent = zrDotNext;
            }

            return new CGStatistics
            {
                AlgorithmName = "Preconditioned Conjugate Gradient",
                HasConverged = false, NumIterationsRequired = maxIterations, NormRatio = resNormRatio
            };
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
