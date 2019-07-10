using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Projection;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;

//TODO: probably needs a builder
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.InterfaceProblem
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
        private readonly IMaxIterationsProvider maxIterationsProvider;
        private readonly IFetiPcgConvergenceFactory pcgConvergenceStrategyFactory;
        private readonly double pcpgConvergenceTolerance;

        private Feti1UnprojectedInterfaceProblemSolver(IMaxIterationsProvider maxIterationsProvider,
            double pcgConvergenceTolerance, IFetiPcgConvergenceFactory pcgConvergenceStrategyFactory)
        {
            this.maxIterationsProvider = maxIterationsProvider;
            this.pcpgConvergenceTolerance = pcgConvergenceTolerance;
            this.pcgConvergenceStrategyFactory = pcgConvergenceStrategyFactory;
        }

        public Vector CalcLagrangeMultipliers(Feti1FlexibilityMatrix flexibility, IFetiPreconditioner preconditioner, 
            Feti1Projection projection, Vector disconnectedDisplacements, Vector rigidBodyModesWork, double globalForcesNorm,
            SolverLogger logger)
        {
            // PCPG starts from the particular lagrange multipliers: λ0 = Q * G * inv(G^T * Q * G) * e
            Vector lagranges = projection.CalcParticularLagrangeMultipliers(rigidBodyModesWork);
            IFetiPcgConvergence pcpgConvergenceStrategy = 
                pcgConvergenceStrategyFactory.CreateConvergenceStrategy(globalForcesNorm);
            var pcpg = new PcpgAlgorithm(maxIterationsProvider, pcpgConvergenceTolerance, pcpgConvergenceStrategy);
            IterativeStatistics stats = 
                pcpg.Solve(flexibility, preconditioner, projection, disconnectedDisplacements, lagranges);
            if (!stats.HasConverged)
            {
                throw new IterativeSolverNotConvergedException(Feti1Solver.name + " did not converge to a solution. PCPG"
                    + $" algorithm run for {stats.NumIterationsRequired} iterations and the residual norm ratio was"
                    + $" {stats.ResidualNormRatioEstimation}");
            }

            // Log statistics about PCG execution
            logger.LogIterativeAlgorithm(stats.NumIterationsRequired, stats.ResidualNormRatioEstimation);
            return lagranges;
        }

        public class Builder
        {
            public IMaxIterationsProvider MaxIterationsProvider { get; set; } = new PercentageMaxIterationsProvider(1.0);
            public IFetiPcgConvergenceFactory PcgConvergenceStrategyFactory { get; set; } =
                new ApproximateResidualConvergence.Factory();
            public double PcgConvergenceTolerance { get; set; } = 1E-7;

            public Feti1UnprojectedInterfaceProblemSolver Build() => new Feti1UnprojectedInterfaceProblemSolver(
                MaxIterationsProvider, PcgConvergenceTolerance, PcgConvergenceStrategyFactory);
        }
    }
}
