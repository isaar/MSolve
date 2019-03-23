using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;

//TODO: probably needs a builder
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti.Feti1
{
    /// <summary>
    /// Uses Preconditioned Conjugate Projected Gradient to solve the original interface linear system. The projection happens
    /// only during PCPG. 
    /// WARNING: PCPG does not converge if the Q matrix used in the projection is not identity. Instead use 
    /// <see cref="Feti1ProjectedInterfaceProblemSolver"/> for heterogeneous or 4th order problems.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Feti1UnprojectedInterfaceProblemSolver : IFeti1InterfaceProblemSolver
    {
        private readonly double pcgConvergenceTolerance;
        private readonly double pcgMaxIterationsOverSize;

        public Feti1UnprojectedInterfaceProblemSolver(double pcgConvergenceTolerance, double pcgMaxIterationsOverSize)
        {
            this.pcgConvergenceTolerance = pcgConvergenceTolerance;
            this.pcgMaxIterationsOverSize = pcgMaxIterationsOverSize;
        }

        public Vector CalcLagrangeMultipliers(Feti1FlexibilityMatrix flexibility, IFetiPreconditioner preconditioner, 
            Feti1Projection projection, Vector disconnectedDisplacements, Vector rigidBodyModesWork, double globalForcesNorm,
            FetiLogger logger)
        {
            // PCPG starts from the particular lagrange multipliers: λ0 = Q * G * inv(G^T * Q * G) * e
            Vector lagranges = projection.CalcParticularLagrangeMultipliers(rigidBodyModesWork);
            var pcpg = new PcpgAlgorithm(pcgMaxIterationsOverSize, pcgConvergenceTolerance, null);
            PcpgStatistics stats = pcpg.Solve(flexibility, preconditioner, projection, disconnectedDisplacements,
                globalForcesNorm, lagranges);

            // Log statistics about PCPG execution
            if (!stats.HasConverged)
            {
                throw new IterativeSolverNotConvergedException(Feti1Solver.name + " did not converge to a solution. PCPG"
                    + $" algorithm run for {stats.NumIterations} iterations and the residual norm ratio was"
                    + $" {stats.ResidualNormEstimateRatio}");
            }
            logger.PcgIterations = stats.NumIterations;
            logger.PcgResidualNormRatio = stats.ResidualNormEstimateRatio;

            return lagranges;
        }
    }
}
