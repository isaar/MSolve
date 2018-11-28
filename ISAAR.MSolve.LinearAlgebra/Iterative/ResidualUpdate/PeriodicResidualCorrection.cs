using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Implement this: "If the tolerance is large, the residual need not be corrected at all" 
//TODO: Add alternative ways to calculate the frequency of corrections.
namespace ISAAR.MSolve.LinearAlgebra.Iterative.ResidualUpdate
{
    /// <summary>
    /// The exact residual (r = b - A*x) is calculated with a fixed frequency to remove the floating point error accumulated by
    /// the more efficient formula used by the iterative algorithm. This approach is presented in section B1 of
    /// "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain", Jonathan Richard Shewchuk, 1994
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PeriodicResidualCorrection: IResidualCorrection
    {
        private IMatrixView matrix;
        private int numIterationsBeforeCorrection; 
        private IVectorView rhs;

        /// <summary>
        /// See <see cref="IResidualCorrection.Initialize(IMatrixView, IVectorView)"/>.
        /// </summary>
        public void Initialize(IMatrixView matrix, IVectorView rhs)
        {
            this.matrix = matrix;
            this.rhs = rhs;
            this.numIterationsBeforeCorrection = (int)Math.Floor(Math.Sqrt(rhs.Length));
        }

        /// <summary>
        /// See <see cref="IResidualCorrection.UpdateResidual(int, IVectorView, IVector, Action{IVector})"/>.
        /// </summary>
        public bool UpdateResidual(int iteration, IVectorView solution, IVector residual, Action<IVector> residualCalculation)
        {
            if ((iteration % numIterationsBeforeCorrection == 0) && (iteration != 0)) //The first iteration uses the correct residual.
            {
                // Calculate the exact residual: r = b - A * x
                residual.CopyFrom(rhs);
                residual.SubtractIntoThis(matrix.MultiplyRight(solution));
                //residual = rhs.Subtract(matrix.MultiplyRight(solution)); //This allocates a new vector r, copies b and GCs the existing r.

                return true;
            }
            else
            {
                // Perform the operation that the iterative algorithm would, without any correction
                residualCalculation(residual);
                return false;
            }
        }
    }
}
