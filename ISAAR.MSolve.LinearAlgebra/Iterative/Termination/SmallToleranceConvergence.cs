using System;
using ISAAR.MSolve.LinearAlgebra.Iterative.ResidualUpdate;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: once the problem with small tolerance has occured once, can we assume that it will always happen? Then we could avoid
//      the regular calculation of the residual and only use r = b - A*x.
namespace ISAAR.MSolve.LinearAlgebra.Iterative.Termination
{
    /// <summary>
    /// Checks if the residual vector has reached the required tolerance. In addition corrects the residual vector, in order to
    /// avoid premature convergence, if the tolerance is close to the limits of the machine's floating point precision, as
    /// described in section B1 of "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain", 
    /// Jonathan Richard Shewchuk, 1994.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SmallToleranceConvergence: IResidualConvergence
    {
        private double limitDotResidual;
        private ILinearTransformation matrix;
        private IVectorView rhs;

        /// <summary>
        /// See 
        /// <see cref="IResidualConvergence.HasConverged(IVectorView, IVector, ref double, bool, Func{IVectorView, double})"/>.
        /// </summary>
        public bool HasConverged(IVectorView solutionVector, IVector residualVector, ref double residualDotProduct,
            bool isResidualCorrected, Func<IVectorView, double> residualDotCalculation)
        {
            // If true, it means convergence, but it may be premature
            bool withinTolerance = residualDotProduct <= limitDotResidual; 

            if (withinTolerance && (!isResidualCorrected)) // The residual may already have been corrected
            {
                // Exact residual: r = b - A * x
                ExactResidual.Calculate(matrix, rhs, solutionVector, residualVector);

                // Recalculate the dot product of the residual vector, as the iterative algorithm would.
                residualDotProduct = residualDotCalculation(residualVector);
                withinTolerance = residualDotProduct <= limitDotResidual;
            }

            return withinTolerance;
        }

        /// <summary>
        /// See <see cref="IResidualConvergence.Initialize(ILinearTransformation, IVectorView, double, double)"/>.
        /// </summary>
        public void Initialize(ILinearTransformation matrix, IVectorView rhsVector, double residualTolerance,
            double initialResidualDotProduct)
        {
            this.matrix = matrix;
            this.rhs = rhsVector;

            // ε^2 * δ0
            // This is more efficient than normalizing and computing square roots. However the order is important, since 
            // ε^2 could be very small and risk precision loss
            this.limitDotResidual = residualTolerance * (residualTolerance * initialResidualDotProduct);
        }
    }
}
