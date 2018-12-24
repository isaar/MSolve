using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;
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
    public class SuiteSparseSolver : ISolver_v2, IDisposable
    {
        private const string name = "SkylineSolver"; // For error messages
        private const bool useSuperNodalFactorization = true; // For faster back/forward substitutions.

        private readonly SymmetricCscAssembler assembler = new SymmetricCscAssembler();
        private readonly IStructuralModel_v2 model;
        private readonly ISubdomain_v2 subdomain;
        private readonly double factorizationPivotTolerance;
        private readonly SuiteSparseSystem linearSystem;

        private bool mustFactorize = true;
        private CholeskySuiteSparse factorization;

        private SuiteSparseSolver(IStructuralModel_v2 model, double factorizationPivotTolerance, IDofOrderer dofOrderer)
        {
            if (model.Subdomains.Count != 1) throw new InvalidSolverException(
                $"{name} can be used if there is only 1 subdomain");
            this.model = model;
            subdomain = model.Subdomains[0];

            linearSystem = new SuiteSparseSystem(subdomain);
            LinearSystems = new Dictionary<int, ILinearSystem_v2>() { { subdomain.ID, linearSystem } };
            linearSystem.MatrixObservers.Add(this);

            this.factorizationPivotTolerance = factorizationPivotTolerance;
            this.DofOrderer = dofOrderer;
        }

        ~SuiteSparseSolver()
        {
            ReleaseResources();
        }

        public IDofOrderer DofOrderer { get; }

        public IReadOnlyDictionary<int, ILinearSystem_v2> LinearSystems { get; }

        public IMatrix BuildGlobalMatrix(ISubdomain_v2 subdomain, IElementMatrixProvider elementMatrixProvider)
            => assembler.BuildGlobalMatrix(subdomain.DofOrdering, subdomain.Elements, elementMatrixProvider);

        /// <summary>
        /// See <see cref="IDisposable.Dispose"/>.
        /// </summary>
        public void Dispose()
        {
            ReleaseResources();
            GC.SuppressFinalize(this);
        }

        public void Initialize() { }

        public void OnMatrixSetting()
        {
            mustFactorize = true;
            if (factorization != null)
            {
                factorization.Dispose();
                factorization = null;
            }
            //TODO: make sure the native memory allocated has been cleared. We need all the available memory we can get.
        }

        /// <summary>
        /// Solves the linear system with back-forward substitution. If the matrix has been modified, it will be refactorized.
        /// </summary>
        public void Solve()
        {
            if (linearSystem.Solution == null) linearSystem.Solution = linearSystem.CreateZeroVector();
            else if (HaveSubdomainDofsChanged()) linearSystem.Solution = linearSystem.CreateZeroVector();
            //else linearSystem.Solution.Clear(); // no need to waste computational time on this

            if (mustFactorize)
            {
                factorization = CholeskySuiteSparse.Factorize(linearSystem.Matrix, useSuperNodalFactorization);
                mustFactorize = false;
            }

            factorization.SolveLinearSystem(linearSystem.RhsVector, linearSystem.Solution);
        }

        //TODO: Create a method in Subdomain (or its DofOrderer) that exposes whether the dofs have changed.
        /// <summary>
        /// The number of dofs might have been changed since the previous Solution vector had been created.
        /// </summary>
        private bool HaveSubdomainDofsChanged() => subdomain.DofOrdering.NumFreeDofs == linearSystem.Solution.Length;

        private void ReleaseResources()
        {
            if (factorization != null)
            {
                factorization.Dispose();
                factorization = null;
            }
        }

        public class Builder
        {
            public Builder() { }

            public IDofOrderer DofOrderer { get; set; }
                = new DofOrderer(new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithSuiteSparseAmd());

            public double FactorizationPivotTolerance { get; set; } = 1E-15;

            public SuiteSparseSolver BuildSolver(IStructuralModel_v2 model)
                => new SuiteSparseSolver(model, FactorizationPivotTolerance, DofOrderer);
        }

        private class SuiteSparseSystem : LinearSystem_v2<SymmetricCscMatrix, Vector>
        {
            internal SuiteSparseSystem(ISubdomain_v2 subdomain) : base(subdomain) { }
            public override Vector CreateZeroVector() => Vector.CreateZero(Subdomain.DofOrdering.NumFreeDofs);
            public override void GetRhsFromSubdomain() => RhsVector = Subdomain.Forces;
        }
    }
}
