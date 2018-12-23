using System;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: perhaps I should name the Action<>
//TODO: perhaps instead of having Initialize() to set the internal state, I should use a factory interface that creates
//      immutable IResidualCorrection objects.
namespace ISAAR.MSolve.LinearAlgebra.Iterative.ResidualUpdate
{
    /// <summary>
    /// Updates the residual vector using an operation provided by the iterative algorithm that uses the
    /// <see cref="IResidualCorrection"/> instance. If necessary, the residual vector is corrected by using the exact formula 
    /// r = b - A*x, instead of the iterative algorithm's. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IResidualCorrection
    {
        /// <summary>
        /// Initializes the internal state of this <see cref="IResidualCorrection"/> instance. Has to be called before each 
        /// linear system solution.
        /// </summary>
        /// <param name="matrix">The linear system's matrix.</param>
        /// <param name="rhs">The linear system's right hand side vector.</param>
        void Initialize(ILinearTransformation matrix, IVectorView rhs);

        /// <summary>
        /// Updates the residual vector and returns true if a correction was necessary. Correcting the residual vector means
        /// using the exact formula r = b - A*x, instead of the iterative algorithm's. 
        /// </summary>
        /// <param name="iteration">The current iteration of the iterative algorithm, starting from 0.</param>
        /// <param name="solution">
        /// The current approximation to the linear system's solution vector. It will not be modified.
        /// </param>
        /// <param name="residual">
        /// The current residual vector. It will be modified as described by <paramref name="residualCalculation"/> or by using
        /// the exact formula r = b - A*x.
        /// </param>
        /// <param name="residualCalculation">
        /// The function used by the iterative algorithm to update the residual vector.
        /// </param>
        bool UpdateResidual(int iteration, IVectorView solution, IVector residual, Action<IVector> residualCalculation);
    }
}
