using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
    /// <summary>
    /// Updates the residual vector r using either the standard operation r = r - α * A*d or the 
    /// exact formula r = b - A*x. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IPcgResidualUpdater
    {
        /// <summary>
        /// Update the residual vector r.
        /// </summary>
        ///<param name="pcg">The Preconditioned Conjugate Gradient algorithm that uses this object.</param>
        /// <param name="residual">The current residual vector r to modify.</param>
        void UpdateResidual(PcgAlgorithmBase pcg, IVector residual);
    }
}
