using System;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: This comment is from when corrections were also applied during checking convergence.This strategy can modify(correct) 
//      r, r*r (CG), s*r (PCG). However that s was computed by preconditioning the uncorrected r. If r is corrected, 
//      should s also be corrected? This is costly, but s will be used extensively afterwards.
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
    public interface IPcgConvergenceStrategy
    {
        /// <summary>
        /// Checks if the residual vector has reached the required tolerance. It may also correct the residual vector and its
        /// dot product with the vector specified by the iterative algorithm.
        /// </summary>
        /// <param name="solutionVector">
        /// The current approximation to the linear system's solution vector x.
        /// </param>
        /// <param name="preconditionedResidualVector">
        /// The current preconditioned residual vector s = inv(M) * r
        /// </param>
        /// <param name="residualDotProduct">
        /// The dot product of the residual vector with the preconditioned residual r^T*inv(M)*r. It is already calculated during 
        /// PCG.
        /// </param>
        bool HasConverged(IVectorView solutionVector, IVectorView preconditionedResidualVector, double residualDotProduct);

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
            double initialResidualDotProduct); //TODO: perhaps the residual tolerance should be injected into the constructor.
    }
}
