using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;

//TODO: Should the matrix managers be injected into the constructor?
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    /// <summary>
    /// The interface problem is solved using PCG. The matrix of the coarse problem KccStar, namely the static condensation of 
    /// the remainder dofs onto the corner dofs is performed explicitly.
    /// </summary>
    public class FetiDPInterfaceProblemSolver : IFetiDPInterfaceProblemSolver
    {
        private readonly IMaxIterationsProvider maxIterationsProvider;
        private readonly double pcgConvergenceTolerance;
        private readonly IFetiPcgConvergenceFactory pcgConvergenceStrategyFactory;

        public FetiDPInterfaceProblemSolver(IMaxIterationsProvider maxIterationsProvider,
            double pcgConvergenceTolerance, IFetiPcgConvergenceFactory pcgConvergenceStrategyFactory)
        {
            this.maxIterationsProvider = maxIterationsProvider;
            this.pcgConvergenceTolerance = pcgConvergenceTolerance;
            this.pcgConvergenceStrategyFactory = pcgConvergenceStrategyFactory;
        }

        public (Vector lagrangeMultipliers, Vector cornerDisplacements) SolveInterfaceProblem(FetiDPFlexibilityMatrix flexibility, 
            IFetiPreconditioner preconditioner, IFetiDPCoarseProblemSolver coarseProblemSolver, 
            Vector globalFcStar, Vector dr, double globalForcesNorm, SolverLogger logger)
        {
            int systemOrder = flexibility.Order;

            // Matrix, preconditioner & rhs
            var pcgMatrix = new InterfaceProblemMatrix(flexibility, coarseProblemSolver);
            var pcgPreconditioner = new InterfaceProblemPreconditioner(preconditioner);
            Vector pcgRhs = CreateInterfaceProblemRhs(flexibility, coarseProblemSolver, globalFcStar, dr);

            // Solve the interface problem using PCG algorithm
            var pcgBuilder = new PcgAlgorithm.Builder();
            pcgBuilder.MaxIterationsProvider = maxIterationsProvider;
            pcgBuilder.ResidualTolerance = pcgConvergenceTolerance;
            pcgBuilder.Convergence = pcgConvergenceStrategyFactory.CreateConvergenceStrategy(globalForcesNorm);
            PcgAlgorithm pcg = pcgBuilder.Build(); //TODO: perhaps use the pcg from the previous analysis if it has reorthogonalization.
            var lagranges = Vector.CreateZero(systemOrder);
            IterativeStatistics stats = pcg.Solve(pcgMatrix, pcgPreconditioner, pcgRhs, lagranges, true,
                () => Vector.CreateZero(systemOrder));

            // Log statistics about PCG execution
            if (!stats.HasConverged)
            {
                throw new IterativeSolverNotConvergedException(FetiDPSolver.name + " did not converge to a solution. PCG"
                    + $" algorithm run for {stats.NumIterationsRequired} iterations and the residual norm ratio was"
                    + $" {stats.ResidualNormRatioEstimation}");
            }
            logger.LogIterativeAlgorithm(stats.NumIterationsRequired, stats.ResidualNormRatioEstimation);

            // Calculate corner displacements: uc = inv(KccStar) * (fcStar + FIrc^T * lagranges)
            Vector uc = flexibility.MultiplyTransposedFIrc(lagranges);
            uc.AddIntoThis(globalFcStar);
            uc = coarseProblemSolver.MultiplyInverseCoarseProblemMatrixTimes(uc);

            return (lagranges, uc);
        }

        private Vector CreateInterfaceProblemRhs(FetiDPFlexibilityMatrix flexibility, 
            IFetiDPCoarseProblemSolver coarseProblemSolver, Vector globalFcStar, Vector dr)
        {
            // rhs = dr - FIrc * inv(KccStar) * fcStar
            Vector rhs = coarseProblemSolver.MultiplyInverseCoarseProblemMatrixTimes(globalFcStar);
            rhs = flexibility.MultiplyFIrc(rhs);
            rhs = dr - rhs;
            return rhs;
        }

        public class Builder
        {
            public IMaxIterationsProvider MaxIterationsProvider { get; set; } = new PercentageMaxIterationsProvider(1.0);
            public IFetiPcgConvergenceFactory PcgConvergenceStrategyFactory { get; set; } =
                new ApproximateResidualConvergence.Factory();
            public double PcgConvergenceTolerance { get; set; } = 1E-7;

            public FetiDPInterfaceProblemSolver Build() => new FetiDPInterfaceProblemSolver(MaxIterationsProvider, 
                PcgConvergenceTolerance, PcgConvergenceStrategyFactory);
        }

        internal class InterfaceProblemMatrix : ILinearTransformation
        {
            private readonly IFetiDPCoarseProblemSolver coarseProblemSolver;
            private readonly FetiDPFlexibilityMatrix flexibility;

            internal InterfaceProblemMatrix(FetiDPFlexibilityMatrix flexibility, IFetiDPCoarseProblemSolver coarseProblemSolver)
            {
                this.flexibility = flexibility;
                this.coarseProblemSolver = coarseProblemSolver;
            }

            public int NumColumns => flexibility.Order;

            public int NumRows => flexibility.Order;

            public void Multiply(IVectorView lhsVector, IVector rhsVector)
            {
                //TODO: remove casts. I think PCG, LinearTransformation and preconditioners should be generic, bounded by 
                //      IVectorView and IVector
                var lhs = (Vector)lhsVector;
                var rhs = (Vector)rhsVector;

                // rhs = (FIrr + FIrc * inv(KccStar) * FIrc^T) * lhs
                rhs.Clear();
                flexibility.MultiplyFIrr(lhs, rhs);
                Vector temp = flexibility.MultiplyTransposedFIrc(lhs);
                temp = coarseProblemSolver.MultiplyInverseCoarseProblemMatrixTimes(temp);
                temp = flexibility.MultiplyFIrc(temp);
                rhs.AddIntoThis(temp);
            }
        }

        internal class InterfaceProblemPreconditioner : IPreconditioner
        {
            private readonly IFetiPreconditioner fetiPreconditioner;

            internal InterfaceProblemPreconditioner(IFetiPreconditioner fetiPreconditioner)
            {
                this.fetiPreconditioner = fetiPreconditioner;
            }

            public void SolveLinearSystem(IVectorView rhsVector, IVector lhsVector)
            {
                //TODO: remove casts. I think PCG, LinearTransformation and preconditioners should be generic, bounded by 
                //      IVectorView and IVector
                var lhs = (Vector)lhsVector;
                var rhs = (Vector)rhsVector;
                fetiPreconditioner.SolveLinearSystem(rhs, lhs);
            }
        }
    }
}
