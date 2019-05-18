using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcpg;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    /// <summary>
    /// The interface problem is solved using PCG. The matrix of the coarse problem KccStar, namely the static condensation of 
    /// the remainder dofs onto the corner dofs is performed explicitly.
    /// </summary>
    public class FetiDPInterfaceProblemSolver : IFetiDPInterfaceProblemSolver
    {
        private readonly PcgAlgorithm.Builder pcgBuilder;
        private CholeskyFull factorizedGlobalKccStar;

        public FetiDPInterfaceProblemSolver(PcgAlgorithm.Builder pcgBuilder)
        {
            this.pcgBuilder = pcgBuilder;
        }

        public void ClearCoarseProblemMatrix()
        {
            factorizedGlobalKccStar = null;
        }

        public void CreateCoarseProblemMatrix(FetiDPDofSeparator dofSeparator, Dictionary<int, CholeskyFull> factorizedKrr, 
            Dictionary<int, Matrix> Krc, Dictionary<int, Matrix> Kcc)
        {
            // Static condensation of remainder dofs (Schur complement).
            var globalKccStar = Matrix.CreateZero(dofSeparator.NumGlobalCornerDofs, dofSeparator.NumGlobalCornerDofs);
            foreach (int s in factorizedKrr.Keys)
            {
                // KccStar[s] = Kcc[s] - Krc[s]^T * inv(Krr[s]) * Krc[s]
                // globalKccStar = sum_over_s(Lc[s]^T * KccStar[s] * Lc[s])
                Matrix Lc = dofSeparator.CornerBooleanMatrices[s];
                Matrix KccStar = Kcc[s] - Krc[s].MultiplyRight(factorizedKrr[s].SolveLinearSystems(Krc[s]), true);
                globalKccStar.AddIntoThis(Lc.ThisTransposeTimesOtherTimesThis(KccStar));
            }

            factorizedGlobalKccStar = globalKccStar.FactorCholesky(true);
        }

        public Vector SolveCoarseProblem(Vector rhs)
        {
            return factorizedGlobalKccStar.SolveLinearSystem(rhs);
        }

        public (Vector lagrangeMultipliers, Vector cornerDisplacements) SolveInterfaceProblem(FetiDPFlexibilityMatrix flexibility, 
            IFetiPreconditioner preconditioner, Vector globalFcStar, Vector dr, double globalForcesNorm, DualSolverLogger logger)
        {
            int systemOrder = flexibility.Order;

            // Matrix and preconditioner
            var pcgMatrix = new InterfaceProblemMatrix(flexibility, factorizedGlobalKccStar);
            var pcgPreconditioner = new InterfaceProblemPreconditioner(preconditioner);

            // rhs = dr - FIrc * inv(KccStar) * fcStar
            Vector pcgRhs = factorizedGlobalKccStar.SolveLinearSystem(globalFcStar);
            pcgRhs = flexibility.MultiplyFIrc(pcgRhs);
            pcgRhs = dr - pcgRhs;

            // Solve the interface problem using PCG algorithm
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
            logger.PcgIterations = stats.NumIterationsRequired;
            logger.PcgResidualNormRatio = stats.ResidualNormRatioEstimation;

            // Calculate corner displacements: uc = inv(KccStar) * (fcStar + FIrc^T * lagranges)
            Vector uc = flexibility.MultiplyTransposedFIrc(lagranges);
            uc.AddIntoThis(globalFcStar);
            uc = factorizedGlobalKccStar.SolveLinearSystem(uc);

            return (lagranges, uc);
        }

        private class InterfaceProblemMatrix : ILinearTransformation
        {
            private readonly FetiDPFlexibilityMatrix flexibility;
            private readonly CholeskyFull factorizedGlobalKccStar;

            internal InterfaceProblemMatrix(FetiDPFlexibilityMatrix flexibility, CholeskyFull factorizedGlobalKccStar)
            {
                this.flexibility = flexibility;
                this.factorizedGlobalKccStar = factorizedGlobalKccStar;
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
                temp = factorizedGlobalKccStar.SolveLinearSystem(temp);
                temp = flexibility.MultiplyFIrc(temp);
                rhs.AddIntoThis(temp);
            }
        }

        private class InterfaceProblemPreconditioner : IPreconditioner
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
