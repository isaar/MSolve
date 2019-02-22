using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.Direct
{
    /// <summary>
    /// Direct solver for models with only 1 subdomain. Uses Cholesky factorization on sparse symmetric positive definite 
    /// matrices stored in symmetric (only the upper triangle) format. Uses native dlls from the SuiteSparse library. 
    /// The factorized matrix and other data are stored in unmanaged memory and properly disposed by this class.
    /// The default behaviour is to apply AMD reordering before Cholesky factorization.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SuiteSparseSolver : SingleSubdomainSolverBase<SymmetricCscMatrix>, IDisposable
    {
        private const bool useSuperNodalFactorization = true; // For faster back/forward substitutions.
        private readonly double factorizationPivotTolerance;

        private bool mustFactorize = true;
        private CholeskySuiteSparse factorization;

        private SuiteSparseSolver(IStructuralModel_v2 model, double factorizationPivotTolerance, IDofOrderer dofOrderer):
            base(model, dofOrderer, new SymmetricCscAssembler(), "SkylineSolver")
        {
            this.factorizationPivotTolerance = factorizationPivotTolerance;
        }

        ~SuiteSparseSolver()
        {
            ReleaseResources();
        }

        public void Dispose()
        {
            ReleaseResources();
            GC.SuppressFinalize(this);
        }

        public override void HandleMatrixWillBeSet()
        {
            mustFactorize = true;
            if (factorization != null)
            {
                factorization.Dispose();
                factorization = null;
            }
            //TODO: make sure the native memory allocated has been cleared. We need all the available memory we can get.
        }

        public override void Initialize() { }

        public override void PreventFromOverwrittingSystemMatrices()
        {
            // The factorization is done over different memory.
        }

        /// <summary>
        /// Solves the linear system with back-forward substitution. If the matrix has been modified, it will be refactorized.
        /// </summary>
        public override void Solve()
        {
            if (linearSystem.Solution == null) linearSystem.Solution = linearSystem.CreateZeroVector();
            //else linearSystem.Solution.Clear(); // no need to waste computational time on this in a direct solver

            if (mustFactorize)
            {
                factorization = CholeskySuiteSparse.Factorize(linearSystem.Matrix, useSuperNodalFactorization);
                mustFactorize = false;
            }

            factorization.SolveLinearSystem(linearSystem.RhsVector, linearSystem.Solution);
        }

        protected override Matrix InverseSystemMatrixTimesOtherMatrix(IMatrixView otherMatrix)
        {
            // Factorization
            if (mustFactorize)
            {
                factorization = CholeskySuiteSparse.Factorize(linearSystem.Matrix, useSuperNodalFactorization);
                mustFactorize = false;
            }

            // Rhs vectors
            Matrix rhsVectors;
            if (otherMatrix is Matrix dense) rhsVectors = dense;
            else rhsVectors = Matrix.CreateFromMatrix(otherMatrix);

            // Solve the linear systems
            Matrix solutionVectors = factorization.SolveLinearSystems(rhsVectors);
            return solutionVectors;
        }

        private void ReleaseResources()
        {
            if (factorization != null)
            {
                factorization.Dispose();
                factorization = null;
            }
        }

        public class Builder : ISolverBuilder
        {
            public Builder() { }

            public IDofOrderer DofOrderer { get; set; }
                = new DofOrderer(new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithSuiteSparseAmd());

            public double FactorizationPivotTolerance { get; set; } = 1E-15;

            ISolver_v2 ISolverBuilder.BuildSolver(IStructuralModel_v2 model) => BuildSolver(model);

            public SuiteSparseSolver BuildSolver(IStructuralModel_v2 model)
                => new SuiteSparseSolver(model, FactorizationPivotTolerance, DofOrderer);
        }
    }
}
