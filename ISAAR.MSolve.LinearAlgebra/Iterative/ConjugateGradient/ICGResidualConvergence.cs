using System;
using System.Collections.Generic;
using System.Text;

//TODO: perhaps instead of having Initialize() to set the internal state, I should use a factory interface that creates
//      immutable IResidualCorrection objects.
//TODO: Join this with the version in PCG
//TODO: ICGConvergence.Initialize() must be called after initializing all properties and before starting the iterations that 
//      will overwrite them. I am not fond of these dependencies.
namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    /// <summary>
    /// Calculates the ratio norm2(f(r)) / norm2(g(r0)), where f and g are vector functions of the residual vector. This ratio
    /// will be used by Conjugate Gradient to check convergence.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface ICGResidualConvergence
    {
        /// <summary>
        /// Calculates the ratio norm2(f(r)) / norm2(g(r0)), where f and g are vector functions of the residual vector.
        /// </summary>
        /// <param name="cg">The Conjugate Gradient Aglorithm that uses this object.</param>
        double EstimateResidualNormRatio(CGAlgorithm cg);

        /// <summary>
        /// Initializes the internal state of this <see cref="ICGResidualConvergence"/> instance. Has to be called immediately 
        /// after calculating the initial residual r0 and r0*r0.
        /// </summary>
        /// <param name="cg">The Conjugate Gradient Aglorithm that uses this object.</param>
        void Initialize(CGAlgorithm cg);
    }
}
