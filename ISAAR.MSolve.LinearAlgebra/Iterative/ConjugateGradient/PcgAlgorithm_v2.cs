using System;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Iterative.ResidualUpdate;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Needs Builder pattern
//TODO: perhaps all quantities should be stored as mutable fields, exposed as readonly properties and the various strategies 
//      should read them from a reference of CG/PCG/PCPG, instead of having them injected.
//TODO: In regular CG, there is a check to perevent premature convergence, by correcting the residual. Can this be done for PCG 
//      as well? Would the preconditioned residual be updated as well?
namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    /// <summary>
    /// Implements the untransformed Preconditioned Conjugate Gradient algorithm for solving linear systems with symmetric 
    /// positive definite matrices. This implementation is based on the algorithm presented in section B3 of 
    /// "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain", Jonathan Richard Shewchuk, 1994
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PcgAlgorithm_v2 : PcgAlgorithmBase_v2
    {
        private readonly IPcgBetaParameterCalculation betaCalculation;

        private PcgAlgorithm_v2(double residualTolerance, IMaxIterationsProvider maxIterationsProvider,
            IPcgConvergenceStrategy pcgConvergence, IResidualCorrection residualCorrection, 
            IPcgBetaParameterCalculation betaCalculation) : 
            base(residualTolerance, maxIterationsProvider, pcgConvergence, residualCorrection)
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
            convergence.Initialize(matrix, rhs, residualTolerance, dotPreconditionedResidualNew);

            // Allocate memory for other vectors, which will be reused during each iteration
            IVector matrixTimesDirection = zeroVectorInitializer();

            for (int iteration = 0; iteration < maxIterations; ++iteration)
            {
                // q = A * d
                matrix.Multiply(direction, matrixTimesDirection);

                // α = δnew / (d * q)
                double stepSize = dotPreconditionedResidualNew / direction.DotProduct(matrixTimesDirection);

                // x = x + α * d
                solution.AxpyIntoThis(direction, stepSize);

                // Normally the residual vector is updated as: r = r - α * q. However corrections might need to be applied.
                residualCorrection.UpdateResidual(iteration, solution, residual,
                    (r) => r.AxpyIntoThis(matrixTimesDirection, -stepSize));

                // s = inv(M) * r
                preconditioner.SolveLinearSystem(residual, preconditionedResidual);

                // δold = δnew
                double dotPreconditionedResidualOld = dotPreconditionedResidualNew;

                // δnew = r * s 
                dotPreconditionedResidualNew = residual.DotProduct(preconditionedResidual);

                // At this point we can check if CG has converged and exit, thus avoiding the uneccesary operations that follow.
                bool hasConverged = convergence.HasConverged(solution, preconditionedResidual, dotPreconditionedResidualNew);
                if (hasConverged)
                {
                    return new CGStatistics
                    {
                        AlgorithmName = "Conjugate Gradient",
                        HasConverged = true,
                        NumIterationsRequired = iteration + 1,
                        ResidualNormRatioEstimation = Math.Sqrt(dotPreconditionedResidualNew) / normResidualInitial //TODO: this is not correct, use the data from convergence strategy
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
                ResidualNormRatioEstimation = Math.Sqrt(dotPreconditionedResidualNew) / normResidualInitial //TODO: this is not correct, use the data from convergence strategy
            };
        }

        /// <summary>
        /// Constructs <see cref="PcgAlgorithm"/> instances, allows the user to specify some or all of the required parameters 
        /// and provides defaults for the rest.
        /// Author: Serafeim Bakalakos
        /// </summary>
        public class Builder : PcgBuilderBase_v2
        {
            /// <summary>
            /// Specifies how to calculate the beta parameter of PCG, which is used to update the direction vector. 
            /// </summary>
            public IPcgBetaParameterCalculation BetaCalculation { get; set; } = new FletcherReevesBeta();

            /// <summary>
            /// Creates a new instance of <see cref="PcgAlgorithm"/>.
            /// </summary>
            public PcgAlgorithm_v2 Build()
            {
                return new PcgAlgorithm_v2(ResidualTolerance, MaxIterationsProvider, Convergence, ResidualCorrection, 
                    BetaCalculation);
            }
        }
    }
}
