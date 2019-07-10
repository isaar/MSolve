using System;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: This comment is from when corrections were also applied during checking convergence.This strategy can modify(correct) 
//      r, r*r (CG), s*r (PCG). However that s was computed by preconditioning the uncorrected r. If r is corrected, 
//      should s also be corrected? This is costly, but s will be used extensively afterwards.
//TODO: perhaps instead of having Initialize() to set the internal state, I should use a factory interface that creates
//      immutable IResidualCorrection objects.
//TODO: Join this with the version in CG
//TODO: ICGConvergence.Initialize() must be called after initializing all properties and before starting the iterations that 
//      will overwrite them. I am not fond of these dependencies.
namespace ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
    /// <summary>
    /// Calculates the ratio norm2(f(r)) / norm2(g(r0)), where f and g are vector functions of the residual vector. This ratio
    /// will be used by Preconditioned Conjugate Gradient to check convergence.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IPcgResidualConvergence
    {
        /// <summary>
        /// Calculates the ratio norm2(f(r)) / norm2(g(r0)), where f and g are vector functions of the residual vector.
        /// </summary>
        ///<param name="pcg">The Preconditioned Conjugate Gradient algorithm that uses this object.</param>
        double EstimateResidualNormRatio(PcgAlgorithmBase pcg);

        /// <summary>
        /// Initializes the internal state of this <see cref="IPcgResidualConvergence"/> instance. Has to be called immediately 
        /// after calculating the initial residual r0 and r0*r0.
        /// </summary>
        ///<param name="pcg">The Preconditioned Conjugate Gradient algorithm that uses this object.</param>
        void Initialize(PcgAlgorithmBase pcg);
    }
}
