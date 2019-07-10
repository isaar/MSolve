using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: rearrange the PCPCG algorithm such that unnecessary operations are avoided.
//TODO: rearrange the PCPCG algorithm, rename variables and use CG utility classes such that PCPCG is consistent with the CG and PCG algorithms.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg
{
    /// <summary>
    /// Implementation of the Preconditioned Conjugate Projected Gradient Method used in the FETI method.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class PcpgAlgorithm
    {
        private const string name = "Preconditioned Conjugate Projected Gradient";
        private readonly IFetiPcgConvergence convergence;
        private readonly IMaxIterationsProvider maxIterationsProvider;
        private readonly double residualTolerance;

        internal PcpgAlgorithm(IMaxIterationsProvider maxIterationsProvider, double residualTolerance, 
            IFetiPcgConvergence convergence)
        {
            this.residualTolerance = residualTolerance;
            this.convergence = convergence;
            this.maxIterationsProvider = maxIterationsProvider;
        }

        internal IterativeStatistics Solve(IInterfaceFlexibilityMatrix matrix, IFetiPreconditioner preconditioner,
            IInterfaceProjection projector, Vector rhs, Vector lagrangeMultipliers)
        {
            int n = matrix.Order;
            int maxIterations = maxIterationsProvider.GetMaxIterations(matrix.Order);

            // r0 = d - F * λ0
            var residual = matrix.Multiply(lagrangeMultipliers);
            residual.LinearCombinationIntoThis(-1.0, rhs, 1.0);

            // Other allocations
            var w = Vector.CreateZero(n);
            var y = Vector.CreateZero(n);
            var z = Vector.CreateZero(n);
            var direction = Vector.CreateZero(n);
            var matrixTimesDirection = Vector.CreateZero(n);
            double residualDotProductPrevious = double.NaN;
            double residualNormRatio = double.NaN;
            double beta = double.NaN;

            for (int iter = 1; iter <= maxIterations; ++iter)
            {
                // w(m-1) = P * r(m-1)
                projector.ProjectVector(residual, w, false);

                // z(m-1) = preconditioner * w(m-1)
                preconditioner.SolveLinearSystem(w, z);

                // Check convergence: usually if ||z|| / ||f|| < tolerance
                residualNormRatio = convergence.EstimateResidualNormRatio(lagrangeMultipliers, z);
                if (residualNormRatio <= residualTolerance)
                {
                    //TODO: is it correct to check for convergence here? How many iterations should I return?
                    return new IterativeStatistics
                    {
                        AlgorithmName = name,
                        HasConverged = true,
                        NumIterationsRequired = iter - 1, //TODO: not sure about this.
                        ResidualNormRatioEstimation = residualNormRatio
                    };
                }

                // y(m-1) = P * z(m-1)
                projector.ProjectVector(z, y, false);

                double residualDotProductCurrent = y.DotProduct(w);

                if (iter == 1)
                {
                    // β(1) = 0
                    beta = 0;

                    // p(1) = y0
                    direction.CopyFrom(y);
                }
                else
                {
                    // β(m) = (y(m-1) * w(m-1)) / (y(m-2) * w(m-2))
                    beta = residualDotProductCurrent / residualDotProductPrevious;

                    // p(m) = y(m-1) + β(m) * p(m-1), if m > 1
                    direction.LinearCombinationIntoThis(beta, y, 1.0);
                }
                residualDotProductPrevious = residualDotProductCurrent;

                // γ(m) = (y(m-1) * w(m-1)) / (p(m) * F * p(m))
                matrix.Multiply(direction, matrixTimesDirection);
                double stepSize = (y * w) / (direction * matrixTimesDirection);

                // λ(m) = λ(m-1) + γ(m) * p(m)
                lagrangeMultipliers.AxpyIntoThis(direction, stepSize);

                // r(m) = r(m-1) -γ(m) * F * p(m)
                residual.AxpyIntoThis(matrixTimesDirection, -stepSize);
            }

            return new IterativeStatistics
            {
                AlgorithmName = name,
                HasConverged = false,
                NumIterationsRequired = maxIterations,
                ResidualNormRatioEstimation = residualNormRatio
            };
        }
    }
}
