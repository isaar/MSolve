using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    /// <summary>
    /// Fletcher-Reeves: beta = (rNew * inv(M) * rNew) / (rOld * inv(M) * rOld).
    /// This is the simplest formula to calculate PCG's beta parameter and does not require any extra memory or calculations.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class FletcherReevesBeta : IPcgBetaParameterCalculation
    {
        /// <summary>
        /// See <see cref="IPcgBetaParameterCalculation.CalculateBeta(IVectorView, IVectorView, double, double)"/>.
        /// </summary>
        public double CalculateBeta(IVectorView residualNew, IVectorView preconditionedResidualNew,
            double dotPreconditionedResidualNew, double dotPreconditionedResidualOld) 
            => dotPreconditionedResidualNew / dotPreconditionedResidualOld;

        /// <summary>
        /// See <see cref="IPcgBetaParameterCalculation.Initialize(IVectorView)"/>.
        /// </summary>
        public void Initialize(IVectorView initialResidual) { } // Do nothing
    }
}
