using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Duplication between this and the CG version
namespace ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
    /// <summary>
    /// Updates the residual vector according to the usual CG formula r = r - α * A*d. No corrections are applied.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class RegularPcgResidualUpdater : IPcgResidualUpdater
    {
        /// <summary>
        /// See <see cref="IPcgResidualUpdater.UpdateResidual(PcgAlgorithmBase, IVector)"/>
        /// </summary>
        public void UpdateResidual(PcgAlgorithmBase pcg, IVector residual)
        {
            // Normally the residual vector is updated as: r = r - α * A*d
            residual.AxpyIntoThis(pcg.MatrixTimesDirection, -pcg.StepSize);
        }
    }
}
