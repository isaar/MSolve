using System;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
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
    public class PCG : PcgBase
    {
        private readonly IPcgBetaParameterCalculation betaCalculation = new FletcherReevesBeta();

        /// <summary>
        /// Initializes a new instance of <see cref="PCG"/> with the specified settings.
        /// </summary>
        /// <param name="maxIterations">The maximum number of iterations before the algorithm terminates.</param>
        /// <param name="residualTolerance">
        /// The algorithm will terminate when sqrt(r*inv(M)*r) / sqrt(r0*inv(M)*r0) &lt;= <paramref name="residualTolerance"/>, 
        /// where x is the current solution vector and x0 the initial guess.
        /// </param>
        public PCG(IMaxIterationsProvider maxIterationsProvider, double residualTolerance): 
            base(maxIterationsProvider, residualTolerance)
        {
        }

        protected override CGStatistics SolveInternal(IMatrixView matrix, IPreconditioner preconditioner, IVectorView rhs, 
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
            betaCalculation.Initialize(residual);
            residualCorrection.Initialize(matrix, rhs);
            residualConvergence.Initialize(matrix, rhs, residualTolerance, dotPreconditionedResidualNew);

            //TODO: Also allocate memory once here, instead of allocating/deallocating it continuously.
            // Declare the vector q = A*d.  
            IVector matrixTimesDirection;

            for (int iteration = 0; iteration < maxIterations; ++iteration)
            {
                // q = A * d
                matrixTimesDirection = matrix.MultiplyRight(direction); //TODO: this allocates a new vector and GCs the existing one

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
    }
}
