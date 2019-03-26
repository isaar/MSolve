using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Implement this: "If the tolerance is large, the residual need not be corrected at all" 
//TODO: Add alternative ways to calculate the frequency of corrections.
namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    /// <summary>
    /// The exact residual (r = b - A*x) is calculated with a fixed frequency to remove the floating point error accumulated by
    /// the more efficient formula used by the iterative algorithm. This approach is presented in section B1 of
    /// "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain", Jonathan Richard Shewchuk, 1994
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PeriodicCorrectionCGResidualUpdater : ICGResidualUpdater
    {
        private int numIterationsBeforeCorrection = int.MinValue;

        /// <summary>
        /// See <see cref="ICGResidualUpdater.UpdateResidual(CGAlgorithm, IVector, out double)"/>
        /// </summary>
        public void UpdateResidual(CGAlgorithm cg, IVector residual, out double resDotRes)
        {
            //TODO: perhaps this should be done in an Initialize() method
            if (numIterationsBeforeCorrection == int.MinValue)
            {
                numIterationsBeforeCorrection = (int)Math.Floor(Math.Sqrt(cg.Rhs.Length));
            }

            if ((cg.Iteration % numIterationsBeforeCorrection == 0) && (cg.Iteration != 0)) //The first iteration uses the correct residual.
            {
                // Calculate the exact residual: r = b - A * x
                ExactResidual.Calculate(cg.Matrix, cg.Rhs, cg.Solution, residual);
            }
            else
            {
                // Normally the residual vector is updated as: r = r - α * A*d
                residual.AxpyIntoThis(cg.MatrixTimesDirection, -cg.StepSize);
            }
            resDotRes = residual.DotProduct(residual);
        }
    }
}
