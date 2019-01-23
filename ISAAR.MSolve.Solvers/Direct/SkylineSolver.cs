using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.Direct
{
    /// <summary>
    /// Direct solver for models with only 1 subdomain. Uses Cholesky factorization on sparse symmetric positive definite 
    /// matrices stored in Skyline format.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SkylineSolver : SingleSubdomainSolverBase<SkylineMatrix>
    {
        private readonly double factorizationPivotTolerance;

        private bool factorizeInPlace = true;
        private bool mustFactorize = true;
        private CholeskySkyline factorizedMatrix;

        private SkylineSolver(IStructuralModel_v2 model, double factorizationPivotTolerance, IDofOrderer dofOrderer):
            base(model, dofOrderer, new SkylineAssembler(), "SkylineSolver")
        {
            this.factorizationPivotTolerance = factorizationPivotTolerance;
        }

        public override void Initialize() { }

        public override void HandleMatrixWillBeSet()
        {
            mustFactorize = true;
            factorizedMatrix = null;
        }

        public override void PreventFromOverwrittingSystemMatrices() => factorizeInPlace = false;

        /// <summary>
        /// Solves the linear system with back-forward substitution. If the matrix has been modified, it will be refactorized.
        /// </summary>
        public override void Solve()
        {
            if (linearSystem.Solution == null) linearSystem.Solution = linearSystem.CreateZeroVector();
            //else linearSystem.Solution.Clear(); // no need to waste computational time on this in a direct solver

            if (mustFactorize)
            {
                factorizedMatrix = linearSystem.Matrix.FactorCholesky(factorizeInPlace, factorizationPivotTolerance); 
                mustFactorize = false;
            }

            factorizedMatrix.SolveLinearSystem(linearSystem.RhsVector, linearSystem.Solution);
        }

        public class Builder
        {
            public Builder() { }

            public IDofOrderer DofOrderer { get; set; } 
                = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public double FactorizationPivotTolerance { get; set; } = 1E-15;

            public SkylineSolver BuildSolver(IStructuralModel_v2 model)
            {
                return new SkylineSolver(model, FactorizationPivotTolerance, DofOrderer);
            }
        }
    }
}
