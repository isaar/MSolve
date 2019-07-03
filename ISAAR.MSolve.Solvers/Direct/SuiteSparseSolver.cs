using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
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

        private SuiteSparseSolver(IStructuralModel model, double factorizationPivotTolerance, IDofOrderer dofOrderer):
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
            var watch = new Stopwatch();
            if (linearSystem.SolutionConcrete == null) linearSystem.SolutionConcrete = linearSystem.CreateZeroVectorConcrete();
            //else linearSystem.Solution.Clear(); // no need to waste computational time on this in a direct solver

            // Factorization
            if (mustFactorize)
            {
                watch.Start();
                factorization = CholeskySuiteSparse.Factorize(linearSystem.Matrix, useSuperNodalFactorization);
                watch.Stop();
                Logger.LogTaskDuration("Matrix factorization", watch.ElapsedMilliseconds);
                watch.Reset();
                mustFactorize = false;
            }

            // Substitutions
            watch.Start();
            factorization.SolveLinearSystem(linearSystem.RhsConcrete, linearSystem.SolutionConcrete);
            watch.Stop();
            Logger.LogTaskDuration("Back/forward substitutions", watch.ElapsedMilliseconds);
            Logger.IncrementAnalysisStep();
        }

        protected override Matrix InverseSystemMatrixTimesOtherMatrix(IMatrixView otherMatrix)
        {
            var watch = new Stopwatch();

            // Factorization
            if (mustFactorize)
            {
                watch.Start();
                factorization = CholeskySuiteSparse.Factorize(linearSystem.Matrix, useSuperNodalFactorization);
                watch.Stop();
                Logger.LogTaskDuration("Matrix factorization", watch.ElapsedMilliseconds);
                watch.Reset();
                mustFactorize = false;
            }

            // Substitutions
            watch.Start();
            Matrix solutionVectors;
            if (otherMatrix is Matrix otherDense) return factorization.SolveLinearSystems(otherDense);
            else
            {
                try
                {
                    // If there is enough memory, copy the RHS matrix to a dense one, to speed up computations. 
                    //TODO: must be benchmarked, if it is actually more efficient than solving column by column.
                    Matrix rhsVectors = otherMatrix.CopyToFullMatrix();
                    solutionVectors = factorization.SolveLinearSystems(rhsVectors);
                }
                catch (InsufficientMemoryException) //TODO: what about OutOfMemoryException?
                {
                    // Solution vectors
                    int systemOrder = linearSystem.Matrix.NumColumns;
                    int numRhs = otherMatrix.NumColumns;
                    solutionVectors = Matrix.CreateZero(systemOrder, numRhs);
                    Vector solutionVector = linearSystem.CreateZeroVectorConcrete();

                    // Solve each linear system separately, to avoid copying the RHS matrix to a dense one.
                    for (int j = 0; j < numRhs; ++j)
                    {
                        if (j != 0) solutionVector.Clear();
                        Vector rhsVector = otherMatrix.GetColumn(j);
                        factorization.SolveLinearSystem(rhsVector, solutionVector);
                        solutionVectors.SetSubcolumn(j, solutionVector);
                    }
                }
            }
            watch.Stop();
            Logger.LogTaskDuration("Back/forward substitutions", watch.ElapsedMilliseconds);
            Logger.IncrementAnalysisStep();
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

            ISolver ISolverBuilder.BuildSolver(IStructuralModel model) => BuildSolver(model);

            public SuiteSparseSolver BuildSolver(IStructuralModel model)
                => new SuiteSparseSolver(model, FactorizationPivotTolerance, DofOrderer);
        }
    }
}
