using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Duplication between this and the CG version
namespace ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
    /// <summary>
    /// The exact residual (r = b - A*x) is calculated with a fixed frequency to remove the floating point error accumulated by
    /// the more efficient formula used by the iterative algorithm. This approach is presented in section B1 of
    /// "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain", Jonathan Richard Shewchuk, 1994
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PeriodicCorrectionPcgResidualUpdater : IPcgResidualUpdater
    {
        private int numIterationsBeforeCorrection = int.MinValue;

        /// <summary>
        /// See <see cref="IPcgResidualUpdater.UpdateResidual(PcgAlgorithmBase, IVector)"/>
        /// </summary>
        public void UpdateResidual(PcgAlgorithmBase pcg, IVector residual)
        {
            //TODO: perhaps this should be done in an Initialize() method
            if (numIterationsBeforeCorrection == int.MinValue)
            {
                numIterationsBeforeCorrection = (int)Math.Floor(Math.Sqrt(pcg.Rhs.Length));
            }

            if ((pcg.Iteration % numIterationsBeforeCorrection == 0) && (pcg.Iteration != 0)) //The first iteration uses the correct residual.
            {
                // Calculate the exact residual: r = b - A * x
                ExactResidual.Calculate(pcg.Matrix, pcg.Rhs, pcg.Solution, residual);
            }
            else
            {
                // Normally the residual vector is updated as: r = r - α * A*d
                residual.AxpyIntoThis(pcg.MatrixTimesDirection, -pcg.StepSize);
            }
        }
    }
}
