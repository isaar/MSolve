using System;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Iterative.ResidualUpdate;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Needs Builder pattern
namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    /// <summary>
    /// Implements the untransformed Preconditioned Conjugate Gradient algorithm for solving linear systems with symmetric 
    /// positive definite matrices. This implementation is based on the algorithm presented in section B3 of 
    /// "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain", Jonathan Richard Shewchuk, 1994
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PcgAlgorithm : PcgAlgorithmBase
    {
        private readonly IPcgBetaParameterCalculation betaCalculation;

        private PcgAlgorithm(double residualTolerance, IMaxIterationsProvider maxIterationsProvider,
            IResidualConvergence residualConvergence, IResidualCorrection residualCorrection, 
            IPcgBetaParameterCalculation betaCalculation) : 
            base(residualTolerance, maxIterationsProvider, residualConvergence, residualCorrection)
        {
            this.betaCalculation = betaCalculation;
        }

        protected override CGStatistics SolveInternal(ILinearTransformation matrix, IPreconditioner preconditioner,
            IVectorView rhs, IVector solution, IVector residual, Func<IVector> zeroVectorInitializer)
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
            betaCalculation.Initialize(residual);
            residualCorrection.Initialize(matrix, rhs);
            residualConvergence.Initialize(matrix, rhs, residualTolerance, dotPreconditionedResidualNew);

            // Allocate memory for other vectors, which will be reused during each iteration
            IVector matrixTimesDirection = rhs.CreateZeroVectorWithSameFormat();

            for (int iteration = 0; iteration < maxIterations; ++iteration)
            {
                // q = A * d
                matrix.Multiply(direction, matrixTimesDirection);

                // α = δnew / (d * q)
                double stepSize = dotPreconditionedResidualNew / direction.DotProduct(matrixTimesDirection);

                // x = x + α * d
                solution.AxpyIntoThis(direction, stepSize);

                // Normally the residual vector is updated as: r = r - α * q. However corrections might need to be applied.
                bool isResidualCorrected = residualCorrection.UpdateResidual(iteration, solution, residual,
                    (r) => r.AxpyIntoThis(matrixTimesDirection, -stepSize));

                // s = inv(M) * r
                preconditioner.SolveLinearSystem(residual, preconditionedResidual);

                // δold = δnew
                double dotPreconditionedResidualOld = dotPreconditionedResidualNew;

                // δnew = r * s 
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

                // The default Fletcher-Reeves formula is: β = δnew / δold = (sNew * rNew) / (sOld * rOld)
                // However we could use a different one, e.g. for variable preconditioning Polak-Ribiere is usually better.
                double beta = betaCalculation.CalculateBeta(residual, preconditionedResidual, 
                    dotPreconditionedResidualNew, dotPreconditionedResidualOld);

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

        /// <summary>
        /// Constructs <see cref="PcgAlgorithm"/> instances, allows the user to specify some or all of the required parameters 
        /// and provides defaults for the rest.
        /// Author: Serafeim Bakalakos
        /// </summary>
        public class Builder : PcgBuilderBase
        {
            /// <summary>
            /// Specifies how to calculate the beta parameter of PCG, which is used to update the direction vector. 
            /// </summary>
            public IPcgBetaParameterCalculation BetaCalculation { get; set; } = new FletcherReevesBeta();

            /// <summary>
            /// Creates a new instance of <see cref="PcgAlgorithm"/>.
            /// </summary>
            public PcgAlgorithm Build()
            {
                return new PcgAlgorithm(ResidualTolerance, MaxIterationsProvider, ResidualConvergence, ResidualCorrection, 
                    BetaCalculation);
            }
        }
    }
}
