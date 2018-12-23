using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Iterative.ResidualUpdate;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: needs to throw exceptions or at least report indefinite, nonsymmetric and singular matrices.
//TODO: initialization should be done by a vector factory, instead of new Vector(..)
namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    /// <summary>
    /// Implements the Conjugate Gradient algorithm for solving linear systems with a positive definite matrix.
    /// The implementation is based on the algorithm presented in section B2 of 
    /// "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain", Jonathan Richard Shewchuk, 1994
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CGAlgorithm
    {
        private readonly IMaxIterationsProvider maxIterationsProvider;
        private readonly IResidualConvergence residualConvergence;
        private readonly IResidualCorrection residualCorrection;
        private readonly double residualTolerance;

        private CGAlgorithm(double residualTolerance, IMaxIterationsProvider maxIterationsProvider,
            IResidualConvergence residualConvergence, IResidualCorrection residualCorrection)
        {
            this.residualTolerance = residualTolerance;
            this.maxIterationsProvider = maxIterationsProvider;
            this.residualConvergence = residualConvergence;
            this.residualCorrection = residualCorrection;
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
            => Solve(new ExplicitMatrixTransformation(matrix), rhs, solution, initialGuessIsZero);

        /// <summary>
        /// Solves the linear system A * x = b, where A = <paramref name="matrix"/> and b = <paramref name="rhs"/>.
        /// Initially x = <paramref name="initialGuess"/> and then it converges to the solution.
        /// </summary>
        /// <param name="matrix">
        /// Represents the matrix A of the linear system A * x = b, which must be symmetric positive definite.
        /// </param>
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
        public CGStatistics Solve(ILinearTransformation matrix, IVectorView rhs, IVector solution, bool initialGuessIsZero) //TODO: find a better way to handle the case x0=0
        {
            //TODO: these will also be checked by the matrix vector multiplication.
            Preconditions.CheckMultiplicationDimensions(matrix.NumColumns, solution.Length);
            Preconditions.CheckSystemSolutionDimensions(matrix.NumRows, rhs.Length);

            // r = b - A * x
            IVector residual;
            if (initialGuessIsZero) residual = rhs.Copy();
            else residual = ExactResidual.Calculate(matrix, rhs, solution);

            return SolveInternal(matrix, rhs, solution, residual);
        }

        private CGStatistics SolveInternal(ILinearTransformation matrix, IVectorView rhs, IVector solution, IVector residual)
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

            // Allocate memory for other vectors, which will be reused during each iteration
            IVector matrixTimesDirection = rhs.CreateZeroVectorWithSameFormat();

            for (int iteration = 0; iteration < maxIterations; ++iteration)
            {
                // q = A * d
                matrix.Multiply(direction, matrixTimesDirection);

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
        

        /// <summary>
        /// Constructs <see cref="CGAlgorithm"/> instances, allows the user to specify some or all of the required parameters and 
        /// provides defaults for the rest.
        /// Author: Serafeim Bakalakos
        /// </summary>
        public class Builder
        {
            /// <summary>
            /// Specifies how to calculate the maximum iterations that the CG algorithm will run for.
            /// </summary>
            public IMaxIterationsProvider MaxIterationsProvider { get; set; } = new PercentageMaxIterationsProvider(1.0);

            /// <summary>
            /// Specifies how the CG algorithm will check that convergence has been reached.
            /// </summary>
            public IResidualConvergence ResidualConvergence { get; set; } = new SimpleConvergence();

            /// <summary>
            /// Specifies how often the residual vector will be corrected by an exact (but costly) calculation.
            /// </summary>
            public IResidualCorrection ResidualCorrection { get; set; } = new NoResidualCorrection();

            /// <summary>
            /// The CG algorithm will converge when norm2(r) / norm2(r0) &lt;= <paramref name="ResidualTolerance"/>, 
            /// where r = A*x is the current residual vector and r0 = A*x0 the initial residual vector.
            /// </summary>
            public double ResidualTolerance { get; set; } = 1E-10;

            /// <summary>
            /// Creates a new instance of <see cref="CGAlgorithm"/>.
            /// </summary>
            public CGAlgorithm Build()
                => new CGAlgorithm(ResidualTolerance, MaxIterationsProvider, ResidualConvergence, ResidualCorrection);
        }
    }
}
