using System;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.ResidualUpdate
{
    /// <summary>
    /// Updates the residual vector according to the iterative algorithm that uses the <see cref="NoResidualCorrection>"/> 
    /// instance. No corrections are applied.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class NoResidualCorrection : IResidualCorrection
    {
        /// <summary>
        /// See <see cref="IResidualCorrection.Initialize(ILinearTransformation, IVectorView)"/>.
        /// </summary>
        public void Initialize(ILinearTransformation matrix, IVectorView rhs)
        {
            // do nothing
        }

        /// <summary>
        /// See <see cref="IResidualCorrection.UpdateResidual(int, IVectorView, IVector, Action{IVector})"/>.
        /// </summary>
        public bool UpdateResidual(int iteration, IVectorView solution, IVector residual, Action<IVector> residualCalculation)
        {
            // Perform the operation that the iterative algorithm would without any correction
            residualCalculation(residual);
            return false;
        }
    }
}
