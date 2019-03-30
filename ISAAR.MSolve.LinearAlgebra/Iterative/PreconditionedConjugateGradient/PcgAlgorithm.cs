using System;
using System.Diagnostics;
using ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Needs Builder pattern
//TODO: perhaps all quantities should be stored as mutable fields, exposed as readonly properties and the various strategies 
//      should read them from a reference of CG/PCG/PCPG, instead of having them injected.
//TODO: In regular CG, there is a check to perevent premature convergence, by correcting the residual. Can this be done for PCG 
//      as well? Would the preconditioned residual be updated as well?
namespace ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
    /// <summary>
    /// Implements the untransformed Preconditioned Conjugate Gradient algorithm for solving linear systems with symmetric 
    /// positive definite matrices. This implementation is based on the algorithm presented in section B3 of 
    /// "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain", Jonathan Richard Shewchuk, 1994
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PcgAlgorithm : PcgAlgorithmBase
    {
        private const string name = "Preconditioned Conjugate Gradient";
        private readonly IPcgBetaParameterCalculation betaCalculation;

        private PcgAlgorithm(double residualTolerance, IMaxIterationsProvider maxIterationsProvider,
            IPcgResidualConvergence pcgConvergence, IPcgResidualUpdater residualUpdater, 
            IPcgBetaParameterCalculation betaCalculation) : 
            base(residualTolerance, maxIterationsProvider, pcgConvergence, residualUpdater)
        {
            this.betaCalculation = betaCalculation;
        }

        protected override IterativeStatistics SolveInternal(int maxIterations, Func<IVector> zeroVectorInitializer)
        {
            // In contrast to the source algorithm, we initialize s here. At each iteration it will be overwritten, 
            // thus avoiding allocating & deallocating a new vector.
            precondResidual = zeroVectorInitializer();

            // d = inv(M) * r
            direction = zeroVectorInitializer();
            Preconditioner.SolveLinearSystem(residual, direction);

            // δnew = δ0 = r * d
            resDotPrecondRes = residual.DotProduct(direction);

            // The convergence and beta strategies must be initialized immediately after the first r and r*inv(M)*r are computed.
            convergence.Initialize(this);
            betaCalculation.Initialize(this);

            // This is also used as output
            double residualNormRatio = double.NaN;

            // Allocate memory for other vectors, which will be reused during each iteration
            matrixTimesDirection = zeroVectorInitializer();

            for (iteration = 0; iteration < maxIterations; ++iteration)
            {
                // q = A * d
                Matrix.Multiply(direction, matrixTimesDirection);

                // α = δnew / (d * q)
                stepSize = resDotPrecondRes / direction.DotProduct(matrixTimesDirection);

                // x = x + α * d
                solution.AxpyIntoThis(direction, stepSize);

                // Normally the residual vector is updated as: r = r - α * q. However corrections might need to be applied.
                residualUpdater.UpdateResidual(this, residual);

                // s = inv(M) * r
                Preconditioner.SolveLinearSystem(residual, precondResidual);

                // δold = δnew
                resDotPrecondResOld = resDotPrecondRes;

                // δnew = r * s 
                resDotPrecondRes = residual.DotProduct(precondResidual);

                // At this point we can check if CG has converged and exit, thus avoiding the uneccesary operations that follow.
                residualNormRatio = convergence.EstimateResidualNormRatio(this);
                Debug.WriteLine($"PCG Iteration = {iteration}: residual norm ratio = {residualNormRatio}");
                if (residualNormRatio <= residualTolerance)
                {
                    return new IterativeStatistics
                    {
                        AlgorithmName = name,
                        HasConverged = true,
                        NumIterationsRequired = iteration + 1,
                        ResidualNormRatioEstimation = residualNormRatio
                    };
                }

                // The default Fletcher-Reeves formula is: β = δnew / δold = (sNew * rNew) / (sOld * rOld)
                // However we could use a different one, e.g. for variable preconditioning Polak-Ribiere is usually better.
                paramBeta = betaCalculation.CalculateBeta(this);

                // d = s + β * d
                //TODO: benchmark the two options to find out which is faster
                //direction = preconditionedResidual.Axpy(direction, beta); //This allocates a new vector d, copies r and GCs the existing d.
                direction.LinearCombinationIntoThis(paramBeta, precondResidual, 1.0); //This performs additions instead of copying and needless multiplications.
            }

            // We reached the max iterations before PCG converged
            return new IterativeStatistics
            {
                AlgorithmName = name,
                HasConverged = false,
                NumIterationsRequired = maxIterations,
                ResidualNormRatioEstimation = residualNormRatio
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
                return new PcgAlgorithm(ResidualTolerance, MaxIterationsProvider, Convergence, ResidualUpdater, 
                    BetaCalculation);
            }
        }
    }
}
