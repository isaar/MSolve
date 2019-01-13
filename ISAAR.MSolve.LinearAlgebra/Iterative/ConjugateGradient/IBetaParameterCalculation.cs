using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    /// <summary>
    /// Defines how the beta parameter of PCG, which is used to update the direction vector, will be calculated. There are
    /// numerous alternatives, which are equivalent only for linear PCG with a constant preconditioner.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IPcgBetaParameterCalculation
    {
        /// <summary>
        /// Initializes the internal state of this <see cref="IPcgBetaParameterCalculation"/> instance. Has to be called once 
        /// per linear system solution.
        /// </summary>
        /// <param name="initialResidual">The first residual vector calculated (r0=b-A*x0).</param>
        void Initialize(IVectorView initialResidual);

        /// <summary>
        /// Calculates the beta parameter of PCG β(t+1), which will be used to update the direction vector d(t+1).
        /// </summary>
        /// <param name="residualNew">The residual vector r(t+1) of the new iteration.</param>
        /// <param name="preconditionedResidualNew">The vector s(t+1) = inv(M) * r(t+1) of the new iteration.</param>
        /// <param name="dotPreconditionedResidualNew">The dot product s(t+1) * r(t+1) of the new iteration.</param>
        /// <param name="dotPreconditionedResidualOld">The dot product s(t) * r(t) of the previous iteration.</param>
        double CalculateBeta(IVectorView residualNew, IVectorView preconditionedResidualNew,
            double dotPreconditionedResidualNew, double dotPreconditionedResidualOld);
    }
}
