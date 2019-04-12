using System;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: once the problem with small tolerance has occured once, can we assume that it will always happen? Then we could avoid
//      the regular calculation of the residual and only use r = b - A*x.
//TODO: This strategy can modify(correct) r, r*r (CG), s*r (PCG). However that s was computed by preconditioning the uncorrected
//      r. If r is corrected, should s also be corrected? This is costly, but s will be used extensively afterwards
//TODO: This behaviour should probably be hardcoded inside the CG/PCG algorithm, since it involves a lot of quantities and other 
//      strategies. It should also happen only when a user defined flag has been set. 
namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    /// <summary>
    /// Updates the residual vector r and the dot product r*r by using the standard formula r = r - α * A*d, checks if CG 
    /// converges and, if it does, recalculates r, r*r using the exact formula: r = b - A*x. This way premature convergence can 
    /// be avoided, if the tolerance is close to the limits of the machine's floating point precision, as described in section
    /// B1 of "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain", Jonathan Richard Shewchuk, 1994.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SmallToleranceCGResidualUpdater : ICGResidualUpdater
    {
        private readonly ICGResidualConvergence convergence;
        private double residualTolerance;

        public SmallToleranceCGResidualUpdater(double residualTolerance, ICGResidualConvergence convergence)
        {
            this.convergence = convergence;
            this.residualTolerance = residualTolerance;
        }

        public void UpdateResidual(CGAlgorithm cg, IVector residual, out double resDotRes)
        {
            // Update the residual vector normally: r = r - α * A*d
            residual.AxpyIntoThis(cg.MatrixTimesDirection, -cg.StepSize);
            resDotRes = residual.DotProduct(residual); //TODO: it is weird that this sets resDotRes and cg.ResDotRes

            // Check if the CG will converge. TODO: Remove duplicate comutations: this check will also be done by CG later.
            double residualNormRatio = convergence.EstimateResidualNormRatio(cg); // let's pray that ICGResidualConvergence does not mutate any fields
            bool hasConverged = residualNormRatio <= residualTolerance;

            // Avoid premature convergence by calculating th exact residual.
            if (hasConverged)
            {
                // Exact residual: r = b - A * x
                ExactResidual.Calculate(cg.Matrix, cg.Rhs, cg.Solution, residual);

                // Recalculate the r * r
                resDotRes = residual.DotProduct(residual);
            }
        }
        
    }
}
