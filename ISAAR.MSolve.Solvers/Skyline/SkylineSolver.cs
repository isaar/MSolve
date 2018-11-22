using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;

//TODO: factorization should be done during Solve() to time correctly. Alternatively, an observer should record the durations.
//TODO: investigate if it is possible to avoid casting the matrix provided by the analyzer/assembler into skyline. Perhaps the
//      the matrix could be obtained directly from the assembler and the analyzer could provide delegates with the operations it 
//      wants done on the matrix, instead of doing them itself.
//TODO: try to abstract the subdomain logic from ther analyzers. I doubt it is possible though.
//TODO: directly pass the single linear system instead of a list that must be checked. The same holds for all solvers and 
//      assemblers.
namespace ISAAR.MSolve.Solvers.Skyline
{
    /// <summary>
    /// Direct solver for models with only 1 subdomain. Uses Cholesky factorization on symmetric positive definite matrices 
    /// stored in Skyline format.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SkylineSolver : ISolver_v2
    {
        private const string name = "SkylineSolver"; // for error messages
        private readonly SkylineAssembler assembler = new SkylineAssembler();
        private readonly IStructuralModel_v2 model;
        private readonly ISubdomain_v2 subdomain;
        private readonly IDofOrderer dofOrderer;
        private readonly double factorizationPivotTolerance;
        private readonly SkylineSystem linearSystem;
        private CholeskySkyline factorizedMatrix;

        public SkylineSolver(IStructuralModel_v2 model, double factorizationPivotTolerance, IDofOrderer dofOrderer)
        {
            if (model.Subdomains.Count != 1) throw new InvalidSolverException(
                $"{name} can be used if there is only 1 subdomain");
            this.model = model;
            this.subdomain = model.Subdomains[0];
            this.linearSystem = new SkylineSystem(subdomain);
            this.LinearSystems = new ILinearSystem_v2[] { linearSystem };

            this.factorizationPivotTolerance = factorizationPivotTolerance;
            this.dofOrderer = dofOrderer;
        }

        public IReadOnlyList<ILinearSystem_v2> LinearSystems { get; }

        public IMatrix BuildGlobalMatrix(ISubdomain_v2 subdomain, IElementMatrixProvider elementMatrixProvider)
            => assembler.BuildGlobalMatrix(subdomain.DofOrdering, subdomain.Elements, elementMatrixProvider);

        //public IMatrix BuildGlobalMatrix(ISubdomain subdomain, IElementMatrixProvider elementMatrixProvider)
        //{
        //    return assembler.BuildGlobalMatrix(subdomain, elementMatrixProvider);
        //}

        public void Initialize()
        {
            // TODO: perhaps order dofs here
        }

        public void OrderDofs() { }
        //public void OrderDofs() => subdomain.DofOrdering = dofOrderer.OrderDofs(model, subdomain);

        /// <summary>
        /// Solves the linear system with back-forward substitution. If the matrix has been modified, it will be refactorized.
        /// </summary>
        public void Solve()
        {
            //TODO: resolve this weird dependency: NewmarkAnalyzer_v2.InitializeInternalVectors() needs the initial solution 
            // (which may be != 0, if it comes from a previous analysis) before the Solver has performed the first system 
            // solution. However, to initialize it we need the rhs vector which is created when the user calls 
            // Model.ConnectDataStructures(). It would be better to only access the number of free dofs, but that also would be 
            // available after Model.ConnectDataStructures().

            if (linearSystem.Solution == null) linearSystem.Solution = linearSystem.CreateZeroVector();
            else if (HasSubdomainDofsChanged()) linearSystem.Solution = linearSystem.CreateZeroVector();
            //else linearSystem.Solution.Clear(); // no need to waste computational time on this

            if (linearSystem.IsMatrixModified)
            {
                factorizedMatrix = linearSystem.Matrix.FactorCholesky(true, factorizationPivotTolerance);
                linearSystem.IsMatrixModified = false; //TODO: this is bad, since someone else might see it as unchanged. Better use observers.
                linearSystem.IsMatrixFactorized = true;
            }

            factorizedMatrix.SolveLinearSystem(linearSystem.RhsVector, linearSystem.Solution);
        }

        //TODO: Create a method in Subdomain (or its DofOrderer) that exposes whether the dofs have changed.
        /// <summary>
        /// The number of dofs might have been changed since the previous Solution vector had been created.
        /// </summary>
        private bool HasSubdomainDofsChanged() => subdomain.DofOrdering.NumFreeDofs == linearSystem.Solution.Length;

        public class Builder
        {
            public Builder() { }

            public IDofOrderer DofOrderer { get; set; } = new SimpleDofOrderer();

            public double FactorizationPivotTolerance { get; set; } = 1E-15;

            public SkylineSolver BuildSolver(IStructuralModel_v2 model) 
                => new SkylineSolver(model, FactorizationPivotTolerance, DofOrderer);
        }

        private class SkylineSystem : LinearSystem_v2<SkylineMatrix, Vector>
        {
            internal SkylineSystem(ISubdomain_v2 subdomain) : base(subdomain) { }
            public override Vector CreateZeroVector() => Vector.CreateZero(Subdomain.DofOrdering.NumFreeDofs);
            public override void GetRhsFromSubdomain() => RhsVector = Subdomain.Forces;
        }
    }
}
