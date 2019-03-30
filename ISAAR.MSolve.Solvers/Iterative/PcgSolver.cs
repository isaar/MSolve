using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.Iterative
{
    /// <summary>
    /// Iterative solver for models with only 1 subdomain. Uses the Proconditioned Conjugate Gradient algorithm.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PcgSolver : SingleSubdomainSolverBase<CsrMatrix>
    {
        private readonly PcgAlgorithm pcgAlgorithm;
        private readonly IPreconditionerFactory preconditionerFactory;

        private bool mustUpdatePreconditioner = true;
        private IPreconditioner preconditioner;

        private PcgSolver(IStructuralModel_v2 model, PcgAlgorithm pcgAlgorithm, IPreconditionerFactory preconditionerFactory, 
            IDofOrderer dofOrderer):
            base(model, dofOrderer, new CsrAssembler(true), "PcgSolver")
        {
            this.pcgAlgorithm = pcgAlgorithm;
            this.preconditionerFactory = preconditionerFactory;
        }

        public override void HandleMatrixWillBeSet()
        {
            mustUpdatePreconditioner = true;
            preconditioner = null;
        }

        public override void Initialize() { }

        public override void PreventFromOverwrittingSystemMatrices()
        {
            // No factorization is done.
        }

        /// <summary>
        /// Solves the linear system with PCG method. If the matrix has been modified, a new preconditioner will be computed.
        /// </summary>
        public override void Solve()
        {
            if (linearSystem.Solution == null) linearSystem.Solution = linearSystem.CreateZeroVector();
            else linearSystem.Solution.Clear();

            if (mustUpdatePreconditioner)
            {
                preconditioner = preconditionerFactory.CreatePreconditionerFor(linearSystem.Matrix);
                mustUpdatePreconditioner = false;
            }

            IterativeStatistics stats = pcgAlgorithm.Solve(linearSystem.Matrix, preconditioner, linearSystem.RhsVector,
                linearSystem.Solution, true, () => linearSystem.CreateZeroVector()); //TODO: This way, we don't know that x0=0, which will result in an extra b-A*0
            if (!stats.HasConverged)
            {
                throw new IterativeSolverNotConvergedException(name + " did not converge to a solution. PCG algorithm run for"
                    + $" {stats.NumIterationsRequired} iterations and the residual norm ratio was"
                    + $" {stats.ResidualNormRatioEstimation}");
            }
        }

        protected override Matrix InverseSystemMatrixTimesOtherMatrix(IMatrixView otherMatrix)
        {
            //TODO: Use a reorthogonalizetion approach when solving multiple rhs vectors. It would be even better if the CG
            //      algorithm exposed a method for solving for multiple rhs vectors.

            // Preconditioning
            if (mustUpdatePreconditioner)
            {
                preconditioner = preconditionerFactory.CreatePreconditionerFor(linearSystem.Matrix);
                mustUpdatePreconditioner = false;
            }

            // Solution vectors
            int systemOrder = linearSystem.Matrix.NumColumns;
            int numRhs = otherMatrix.NumColumns;
            var solutionVectors = Matrix.CreateZero(systemOrder, numRhs);
            Vector solutionVector = linearSystem.CreateZeroVector();

            // Solve each linear system
            for (int j = 0; j < numRhs; ++j)
            {
                if (j != 0) solutionVector.Clear();

                //TODO: we should make sure this is the same type as the vectors used by this solver, otherwise vector operations
                //      in CG will be slow.
                Vector rhsVector = otherMatrix.GetColumn(j);

                IterativeStatistics stats = pcgAlgorithm.Solve(linearSystem.Matrix, preconditioner, rhsVector,
                    solutionVector, true, () => linearSystem.CreateZeroVector());

                solutionVectors.SetSubcolumn(j, solutionVector);
            }

            return solutionVectors;
        }

        public class Builder : ISolverBuilder
        {
            public IDofOrderer DofOrderer { get; set; }
                = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public PcgAlgorithm PcgAlgorithm { get; set; } = (new PcgAlgorithm.Builder()).Build();

            public IPreconditionerFactory PreconditionerFactory { get; set; } = new JacobiPreconditioner.Factory();

            ISolver_v2 ISolverBuilder.BuildSolver(IStructuralModel_v2 model) => BuildSolver(model);

            public PcgSolver BuildSolver(IStructuralModel_v2 model) 
                => new PcgSolver(model, PcgAlgorithm, PreconditionerFactory, DofOrderer);
        }
    }
}
