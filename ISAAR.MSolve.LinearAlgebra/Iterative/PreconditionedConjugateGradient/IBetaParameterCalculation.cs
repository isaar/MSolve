using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
    /// <summary>
    /// Defines how the beta parameter of PCG, which is used to update the direction vector, will be calculated. There are
    /// numerous alternatives, which are equivalent only for linear PCG with a constant preconditioner.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IPcgBetaParameterCalculation
    {
        /// <summary>
        /// Calculates the beta parameter of PCG β(t+1), which will be used to update the direction vector d(t+1).
        /// </summary>
        ///<param name="pcg">The Preconditioned Conjugate Gradient algorithm that uses this object.</param>
        double CalculateBeta(PcgAlgorithmBase pcg);

        /// <summary>
        /// Initializes the internal state of this <see cref="IPcgBetaParameterCalculation"/> instance. Has to be called 
        /// immediately after calculating the initial residual r0 and r0*r0.
        /// </summary>
        ///<param name="pcg">The Preconditioned Conjugate Gradient algorithm that uses this object.</param>
        void Initialize(PcgAlgorithmBase pcg);
    }
}
