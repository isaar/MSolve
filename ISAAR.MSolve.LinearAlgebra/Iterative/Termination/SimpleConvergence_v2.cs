using System;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.Termination
{
    /// <summary>
    /// Checks if the residual vector has reached the required tolerance, without performing corrections.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SimpleConvergence_v2 : IPcgConvergenceStrategy
    {
        private double residualDotProductLimit;

        /// <summary>
        /// See 
        /// <see cref="IResidualConvergence.HasConverged(IVectorView, IVector, ref double, bool, Func{IVectorView, double})"/>.
        /// </summary>
        public bool HasConverged(IVectorView solutionVector, IVectorView preconditionedResidualVector, double residualDotProduct)
            => residualDotProduct <= residualDotProductLimit;

        /// <summary>
        /// See <see cref="IResidualConvergence.Initialize(ILinearTransformation, IVectorView, double, double)"/>.
        /// </summary>
        public void Initialize(ILinearTransformation matrix, IVectorView rhsVector, double residualTolerance,
            double initialResidualDotProduct)
        {
            // ε^2 * δ0
            // This is more efficient than normalizing and computing square roots. However the order is important, since 
            // ε^2 could be very small and risk precision loss
            this.residualDotProductLimit = residualTolerance * (residualTolerance * initialResidualDotProduct);
        }
    }
}
