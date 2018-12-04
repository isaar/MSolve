using System;
using System.Diagnostics;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: needs builder
namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    /// <summary>
    /// Implements the untransformed Preconditioned Conjugate Gradient algorithm for solving linear systems with symmetric 
    /// positive definite matrices. The implementation is based on the algorithm presented in pages 51-54 of the PhD dissertation 
    /// "Seismic soil-structure interaction with finite elements and the method of substructures", George Stavroulakis, 2014
    /// Authors: Serafeim Bakalakos, George Stavroulakis 
    /// </summary>
    public class PcgWithReorthogonalization: PcgBase
    {
        //TODO: this could be abstracted to use a cyclic cache.
        private readonly PcgReorthogonalizationCache reorthoCache = new PcgReorthogonalizationCache();

        /// <summary>
        /// Initializes a new instance of <see cref="PcgWithReorthogonalization"/> with the specified settings.
        /// </summary>
        /// <param name="maxIterations">The maximum number of iterations before the algorithm terminates.</param>
        /// <param name="residualTolerance">
        /// The algorithm will terminate when sqrt(r*inv(M)*r) / sqrt(r0*inv(M)*r0) &lt;= <paramref name="residualTolerance"/>, 
        /// where x is the current solution vector and x0 the initial guess.
        /// </param>
        public PcgWithReorthogonalization(IMaxIterationsProvider maxIterationsProvider, double residualTolerance) :
            base(maxIterationsProvider, residualTolerance)
        {
        }

        /// <summary>
        /// Calculates the initial approximation to the linear system's solution vector, by using a series of conjugate direction 
        /// vectors that have been stored previously by PCG during the solution of other linear systems with the same matrix. 
        /// This method should be used to solve linear systems with different right hand side vectors and the same matrix.
        /// </summary>
        /// <param name="rhsNew">The right hand side vector of the new linear system.</param>
        /// <param name="initialSolution">
        /// The initial approximation to the solution vector, which PCG will improve. It will be overwritten by this method.
        /// </param>
        /// <param name="isSolutionZero">
        /// Set to true if <paramref name="initialSolution"/> is the zero vector, to avoid clearing it.
        /// </param>
        /// <exception cref="InvalidOperationException">Thrown if there are no direction vectors stored yet.</exception>
        public void CalculateInitialSolutionFromStoredDirections(IVectorView rhsNew, IVector initialSolution, 
            bool isSolutionZero)
        {
            if (reorthoCache.Directions.Count < 1) throw new InvalidOperationException("There are no direction vectors stored.");
            if (!isSolutionZero) initialSolution.Clear();

            //TODO: An implementation by G. Stavroulakis discarded the last stored direction vector at this point. Why?
            //reorthoCache.RemoveNewDirectionVectorData(1);

            // x0 = D_nd * x_d, x_d = inv(Q_nd * D_nd) * D_nd^T * b
            // D_nd = [d_1 ... d_nd], Q_nd = A * D_nd = [q_1 ... q_nd], Q_nd * D_nd = diag([d1*A*d1 ... d_nd*A*d_nd])
            for (int i = 0; i < reorthoCache.Directions.Count; ++i)
            {
                // x_d[i] = 1/(d_i * q_i) * (d_i * b)
                double xd = reorthoCache.Directions[i].DotProduct(rhsNew) / reorthoCache.DirectionsTimesMatrixTimesDirections[i];

                Debug.Assert(!double.IsNaN(xd));
                Debug.Assert(!double.IsPositiveInfinity(xd));
                Debug.Assert(!double.IsNegativeInfinity(xd));

                // x0 += d_i * x_d[i]
                initialSolution.AxpyIntoThis(reorthoCache.Directions[i], xd);
            }
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

            // q = A * d
            IVector matrixTimesDirection = matrix.MultiplyRight(direction);
            double directionTimesMatrixTimesDirection = direction.DotProduct(matrixTimesDirection);

            // Update the direction vectors cache
            reorthoCache.StoreDirectionData(direction, matrixTimesDirection, directionTimesMatrixTimesDirection);

            // δnew = δ0 = r * d
            double dotPreconditionedResidualNew = residual.DotProduct(direction);

            // This is only used as output
            double normResidualInitial = Math.Sqrt(dotPreconditionedResidualNew);

            // Initialize the strategy objects
            residualCorrection.Initialize(matrix, rhs);
            residualConvergence.Initialize(matrix, rhs, residualTolerance, dotPreconditionedResidualNew);

            //TODO: Find proof that this correct. Why is it better than the default formula α = (r * s) / (d * q)?
            // α = (d * r) / (d * q) = (d * r) / (d * (A * d)) 
            // In the first iteration all multiplications have already been performed.
            double stepSize = dotPreconditionedResidualNew / directionTimesMatrixTimesDirection;

            for (int iteration = 0; iteration < maxIterations; ++iteration)
            {
                // x = x + α * d
                solution.AxpyIntoThis(direction, stepSize);

                // Normally the residual vector is updated as: r = r - α * q. However corrections might need to be applied.
                bool isResidualCorrected = residualCorrection.UpdateResidual(iteration, solution, residual,
                    (r) => r.AxpyIntoThis(matrixTimesDirection, -stepSize));

                // s = inv(M) * r
                preconditioner.SolveLinearSystem(residual, preconditionedResidual);

                // δold = δnew
                double dotPreconditionedResidualOld = dotPreconditionedResidualNew;

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

                // Update the direction vector using previous cached direction vectors.
                UpdateDirectionVector(preconditionedResidual, direction);

                // q = A * d
                matrixTimesDirection = matrix.MultiplyRight(direction); //TODO: this allocates a new vector and GCs the existing one
                directionTimesMatrixTimesDirection = direction.DotProduct(matrixTimesDirection);

                // Update the direction vectors cache
                reorthoCache.StoreDirectionData(direction, matrixTimesDirection, directionTimesMatrixTimesDirection);

                //TODO: Find proof that this correct. Why is it better than the default formula α = (r * s) / (d * q)?
                // α = (d * r) / (d * q) = (d * r) / (d * (A * d)) 
                stepSize = direction.DotProduct(residual) / directionTimesMatrixTimesDirection;
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

        private void UpdateDirectionVector(IVectorView preconditionedResidual, IVector direction)
        {
            // d = s - sum(β_i * d_i), 0 <= i < currentIteration
            // β_i = (s * q_i) / (d_i * q_i)
            direction.CopyFrom(preconditionedResidual);
            for (int i = 0; i < reorthoCache.Directions.Count; ++i)
            {
                double beta = preconditionedResidual.DotProduct(reorthoCache.MatrixTimesDirections[i])
                    / reorthoCache.DirectionsTimesMatrixTimesDirections[i];
                direction.AxpyIntoThis(reorthoCache.Directions[i], -beta);
            }
        }
    }
}
