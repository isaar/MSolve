using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: See B5 of painless CG. We can avoid a vector copy per iteration, by doing the old-new dot product before updating s.
//      Not sure if that is correct.
//TODO: Perhaps the previous residual should be stored here instead of being exposed by the PCG algorithm itself.
namespace ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
    /// <summary>
    /// Polak-Ribiere: beta = max{((rNew-rOld) * inv(M) * rNew) / (rOld * inv(M) * rOld), 0}.
    /// This formula usually improves performance for variable preconditioners (changing between iterations). However it 
    /// requires an extra vector to be copied and stored at each iteration.
    /// </summary>
    public class PolakRibiereBeta : IPcgBetaParameterCalculation
    {
        private IVector residualOld;

        /// <summary>
        /// See <see cref="IPcgBetaParameterCalculation.Initialize(IVectorView)"/>.
        /// </summary>
        public void Initialize(PcgAlgorithmBase pcg) => residualOld = pcg.Residual.Copy();

        /// <summary>
        /// See <see cref="IPcgBetaParameterCalculation.CalculateBeta(PcgAlgorithmBase)"/>.
        /// </summary>
        public double CalculateBeta(PcgAlgorithmBase pcg)
        {
            double nominator = pcg.ResDotPrecondRes - pcg.PrecondResidual.DotProduct(residualOld);
            residualOld.CopyFrom(pcg.Residual);
            return (nominator > 0.0) ? nominator / pcg.ResDotPrecondResOld : 0.0;
        }
    }
}
