using System;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: This strategy can modify(correct) r, r*r (CG), s*r (PCG). However that s was computed by preconditioning the uncorrected
//      r. If r is corrected, should s also be corrected? This is costly, but s will be used extensively afterwards.
//TODO: perhaps I should name the Func<>
//TODO: perhaps instead of having Initialize() to set the internal state, I should use a factory interface that creates
//      immutable IResidualCorrection objects.
namespace ISAAR.MSolve.LinearAlgebra.Iterative.Termination
{
    /// <summary>
    /// Checks if the residual vector has reached the required tolerance. It may also correct the residual vector and its
    /// dot product with the vector specified by the iterative algorithm that uses this <see cref="IResidualConvergence"/>, in 
    /// order to avoid numerical complications. The correct residual vector is calculated as r = b - A * x.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IResidualConvergence
    {
        /// <summary>
        /// Checks if the residual vector has reached the required tolerance. It may also correct the residual vector and its
        /// dot product with the vector specified by the iterative algorithm.
        /// </summary>
        /// <param name="solutionVector">
        /// The current approximation to the linear system's solution vector. It will not be modified.
        /// </param>
        /// <param name="residualVector">
        /// The current residual vector. It may be modified by using the exact formula r = b - A*x.
        /// </param>
        /// <param name="residualDotProduct">
        /// The dot product of the residual vector, with a vector specified by the iterative algorithm. If
        /// <paramref name="residualVector"/> is corrected, <paramref name="residualDotProduct"/> will also be modified according
        /// to <paramref name="residualDotCalculation"/>.
        /// </param>
        /// <param name="isResidualCorrected">
        /// True if the residual vector has already been corrected during the current iteration of the iterative algorithm.
        /// </param>
        /// <param name="residualDotCalculation">
        /// The function used by the iterative algorithm to calculate the dot product of the residual vector with another
        /// vector, specified by the iterative algorithm.
        /// </param>
        bool HasConverged(IVectorView solutionVector, IVector residualVector, ref double residualDotProduct,
            bool isResidualCorrected, Func<IVectorView, double> residualDotCalculation);

        /// <summary>
        /// Initializes the internal state of this <see cref="IResidualConvergence"/> instance. Has to be called before each 
        /// linear system solution.
        /// </summary>
        /// <param name="matrix">The linear system's matrix.</param>
        /// <param name="rhsVector">The linear system's right hand side vector.</param>
        /// <param name="residualTolerance">
        /// The max allowable value of the residual vector's normalized dot product, to consider that the iterative algorithm 
        /// has converged.
        /// </param>
        /// <param name="initialResidualDotProduct">
        /// The initial dot product of the residual vector, with a vector specified by the iterative algorithm.
        /// </param>
        void Initialize(ILinearTransformation matrix, IVectorView rhsVector, double residualTolerance, 
            double initialResidualDotProduct);
    }
}
