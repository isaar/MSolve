using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;

//TODO: Improve CG, PCG with strategy patterns(for seach directions, beta calculation, etc), avoid the first r=b-A*0 
//TODO: IIndexable2D is not a good choice if all solvers must cast it to the matrix types the operate on.
//TODO: perhaps the internal vectors of PCG can be cleared and reused.
namespace ISAAR.MSolve.Solvers.PCG
{
    /// <summary>
    /// Iterative solver for models with only 1 subdomain. Uses the Proconditioned Conjugate Gradient algorithm.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PcgSolver : ISolver_v2
    {
        private const string name = "PcgSolver"; // for error messages
        private readonly CsrAssembler assembler = new CsrAssembler(true);
        private readonly IStructuralModel_v2 model;
        private readonly ISubdomain_v2 subdomain;
        private readonly CsrSystem linearSystem;
        private readonly PcgAlgorithm pcgAlgorithm;
        private readonly IPreconditionerFactory preconditionerFactory;
        private IPreconditioner preconditioner;

        public PcgSolver(IStructuralModel_v2 model, PcgAlgorithm pcgAlgorithm, IPreconditionerFactory preconditionerFactory, 
            IDofOrderer dofOrderer)
        {
            if (model.Subdomains.Count != 1) throw new InvalidSolverException(
                $"{name} can be used if there is only 1 subdomain");
            this.model = model;
            this.subdomain = model.Subdomains[0];
            this.linearSystem = new CsrSystem(subdomain);
            this.LinearSystems = new ILinearSystem_v2[] { linearSystem };

            this.pcgAlgorithm = pcgAlgorithm;
            this.preconditionerFactory = preconditionerFactory;
            this.DofOrderer = dofOrderer;
        }

        public IDofOrderer DofOrderer { get; }

        public IReadOnlyList<ILinearSystem_v2> LinearSystems { get; }

        public IMatrix BuildGlobalMatrix(ISubdomain_v2 subdomain, IElementMatrixProvider elementMatrixProvider)
            => assembler.BuildGlobalMatrix(subdomain.DofOrdering, subdomain.Elements, elementMatrixProvider);

        public void Initialize() { }

        /// <summary>
        /// Solves the linear system with PCG method. If the matrix has been modified, a new preconditioner will be computed.
        /// </summary>
        public void Solve()
        {
            if (linearSystem.Solution == null) linearSystem.Solution = linearSystem.CreateZeroVector();
            else if (HasSubdomainDofsChanged()) linearSystem.Solution = linearSystem.CreateZeroVector();
            else linearSystem.Solution.Clear(); // In iterative algorithms we initialize the solution vector to 0.

            if (linearSystem.IsMatrixModified)
            {
                preconditioner = preconditionerFactory.CreatePreconditionerFor(linearSystem.Matrix);
                linearSystem.IsMatrixModified = false;
            }

            CGStatistics stats = pcgAlgorithm.Solve(linearSystem.Matrix, preconditioner, linearSystem.RhsVector,
                linearSystem.Solution, true, () => linearSystem.CreateZeroVector()); //TODO: This way, we don't know that x0=0, which will result in an extra b-A*0
        }

        //TODO: Create a method in Subdomain (or its DofOrderer) that exposes whether the dofs have changed.
        private bool HasSubdomainDofsChanged() => subdomain.DofOrdering.NumFreeDofs == linearSystem.Solution.Length;

        public class Builder
        {
            private readonly PcgAlgorithm pcg;

            public Builder(PcgAlgorithm pcg) => this.pcg = pcg;

            public IDofOrderer DofOrderer { get; set; } = new SimpleDofOrderer();
        
            public IPreconditionerFactory PreconditionerFactory { get; set; } = new JacobiPreconditioner.Factory();

            public PcgSolver BuildSolver(IStructuralModel_v2 model)
                => new PcgSolver(model, pcg, PreconditionerFactory, DofOrderer);
        }

        private class CsrSystem : LinearSystem_v2<CsrMatrix, Vector>
        {
            internal CsrSystem(ISubdomain_v2 subdomain) : base(subdomain) { }
            public override Vector CreateZeroVector() => Vector.CreateZero(Subdomain.DofOrdering.NumFreeDofs);
            public override void GetRhsFromSubdomain() => RhsVector = Subdomain.Forces;
        }
    }
}
