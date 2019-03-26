using System;
using System.Diagnostics;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: I would rather implement reorthogonalization as an alternative strategy, rather than a different class.
//TODO: needs builder
namespace ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
    /// <summary>
    /// Implements the untransformed Preconditioned Conjugate Gradient algorithm for solving linear systems with symmetric 
    /// positive definite matrices. The implementation is based on the algorithm presented in pages 51-54 of the PhD dissertation 
    /// "Seismic soil-structure interaction with finite elements and the method of substructures", George Stavroulakis, 2014
    /// Authors: Serafeim Bakalakos, George Stavroulakis 
    /// </summary>
    public class PcgWithReorthogonalization : PcgAlgorithmBase
    {
        private const string name = "PCG with reorthogonalization";

        //TODO: this could be abstracted to use a cyclic cache.
        private readonly PcgReorthogonalizationCache reorthoCache = new PcgReorthogonalizationCache();

        private PcgWithReorthogonalization(double residualTolerance, IMaxIterationsProvider maxIterationsProvider,
            IPcgResidualConvergence residualConvergence, IPcgResidualUpdater residualCorrection) :
            base(residualTolerance, maxIterationsProvider, residualConvergence, residualCorrection)
        {
        }

        /// <summary>
        /// The dot product d * (A*d), where d is the direction vector <see cref="PcgAlgorithmBase.Direction"/>.
        /// </summary>
        public double DirectionTimesMatrixTimesDirection { get; private set; }

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

        /// <summary>
        /// See <see cref="PcgAlgorithmBase.Clear"/>
        /// </summary>
        public override void Clear()
        {
            base.Clear();
            DirectionTimesMatrixTimesDirection = 0.0;
        }

        protected override IterativeStatistics SolveInternal(int maxIterations, Func<IVector> zeroVectorInitializer)
        {
            // In contrast to the source algorithm, we initialize s here. At each iteration it will be overwritten, 
            // thus avoiding allocating deallocating a new vector.
            precondResidual = zeroVectorInitializer();

            // d = inv(M) * r
            direction = zeroVectorInitializer();
            Preconditioner.SolveLinearSystem(residual, direction);

            // q = A * d
            matrixTimesDirection = zeroVectorInitializer();
            Matrix.Multiply(direction, matrixTimesDirection);
            double directionTimesMatrixTimesDirection = direction.DotProduct(matrixTimesDirection);

            // Update the direction vectors cache
            reorthoCache.StoreDirectionData(this);

            // δnew = δ0 = r * d
            double resDotPrecondRes = residual.DotProduct(direction);

            // The convergence strategy must be initialized immediately after the first r and r*inv(M)*r are computed.
            convergence.Initialize(this);

            // This is also used as output
            double residualNormRatio = double.NaN;

            //TODO: Find proof that this correct. Why is it better than the default formula α = (r * s) / (d * q)?
            // α = (d * r) / (d * q) = (d * r) / (d * (A * d)) 
            // In the first iteration all multiplications have already been performed.
            stepSize = resDotPrecondRes / directionTimesMatrixTimesDirection;

            for (int iteration = 0; iteration < maxIterations; ++iteration)
            {
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

                /// At this point we can check if CG has converged and exit, thus avoiding the uneccesary operations that follow.
                residualNormRatio = convergence.EstimateResidualNormRatio(this);
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

                // Update the direction vector using previous cached direction vectors.
                UpdateDirectionVector(precondResidual, direction);

                // q = A * d
                Matrix.Multiply(direction, matrixTimesDirection);
                directionTimesMatrixTimesDirection = direction.DotProduct(matrixTimesDirection);

                // Update the direction vectors cache
                reorthoCache.StoreDirectionData(this);

                //TODO: Find proof that this correct. Why is it better than the default formula α = (r * s) / (d * q)?
                // α = (d * r) / (d * q) = (d * r) / (d * (A * d)) 
                stepSize = direction.DotProduct(residual) / directionTimesMatrixTimesDirection;
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

        /// <summary>
        /// Constructs <see cref="PcgWithReorthogonalization"/> instances, allows the user to specify some or all of the 
        /// required parameters and provides defaults for the rest.
        /// Author: Serafeim Bakalakos
        /// </summary>
        public class Builder: PcgBuilderBase
        {
            /// <summary>
            /// Creates a new instance of <see cref="PcgWithReorthogonalization"/>.
            /// </summary>
            public PcgWithReorthogonalization Build()
            {
                return new PcgWithReorthogonalization(ResidualTolerance, MaxIterationsProvider, Convergence, ResidualUpdater);
            }
        }
    }
}
