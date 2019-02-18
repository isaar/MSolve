using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: See B5 of painless CG. We can avoid a vector copy per iteration, by doing the old-new dot product before updating s.
//      Not sure if that is correct.
//TODO: Perhaps some other strategy also needs to store the previous residual. We should store and copy it centrally. 
namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    /// <summary>
    /// Polak-Ribiere: beta = max{((rNew-rOld) * inv(M) * rNew) / (rOld * inv(M) * rOld), 0}.
    /// This formula usually improves performance for variable preconditioners (changing between iterations). However it 
    /// requires an extra vector to be copied and stored at each iteration.
    /// </summary>
    public class PolakRibiereBeta : IPcgBetaParameterCalculation
    {
        private IVectorView residualOld;

        /// <summary>
        /// See <see cref="IPcgBetaParameterCalculation.CalculateBeta(IVectorView, IVectorView, double, double)"/>.
        /// </summary>
        public double CalculateBeta(IVectorView residualNew, IVectorView preconditionedResidualNew, 
            double dotPreconditionedResidualNew, double dotPreconditionedResidualOld)
        {
            double nominator = dotPreconditionedResidualNew - preconditionedResidualNew.DotProduct(residualOld);
            residualOld = residualNew.Copy();
            return (nominator > 0.0) ? nominator / dotPreconditionedResidualOld : 0.0;
        }

        /// <summary>
        /// See <see cref="IPcgBetaParameterCalculation.Initialize(IVectorView)"/>.
        /// </summary>
        public void Initialize(IVectorView initialResidual) => residualOld = initialResidual.Copy();
    }
}
