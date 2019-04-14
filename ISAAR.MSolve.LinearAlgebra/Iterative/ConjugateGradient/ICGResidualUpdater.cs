using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    /// <summary>
    /// Updates the residual vector r and the dot product r*r by using either the standard operation r = r - α * A*d or the 
    /// exact formula r = b - A*x. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface ICGResidualUpdater
    {
        /// <summary>
        /// Update the residual vector r and the dot product r*r
        /// </summary>
        /// <param name="cg">The Conjugate Gradient Aglorithm that uses this object.</param>
        /// <param name="residual">The current residual vector r to modify.</param>
        /// <param name="resDotRes">The current product r*r to modify.</param>
        void UpdateResidual(CGAlgorithm cg, IVector residual, out double resDotRes);
    }
}
