using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.Algorithms.CG;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;

//TODO: Improve CG, PCG with strategy patterns(for seach directions, beta calculation, etc), avoid the first r=b-A*0 
//TODO: IIndexable2D is not a good choice if all solvers must cast it to the matrix types the operate on.
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
        private readonly IDofOrderer dofOrderer;
        private readonly CsrSystem linearSystem;
        private readonly PreconditionedConjugateGradient pcgAlgorithm;
        private readonly IPreconditionerFactory preconditionerFactory;
        private IPreconditioner preconditioner;

        public PcgSolver(IStructuralModel_v2 model, PreconditionedConjugateGradient pcgAlgorithm, 
            IPreconditionerFactory preconditionerFactory, IDofOrderer dofOrderer)
        {
            if (model.Subdomains.Count != 1) throw new InvalidSolverException(
                $"{name} can be used if there is only 1 subdomain");
            this.model = model;
            this.subdomain = model.Subdomains[0];
            this.linearSystem = new CsrSystem(subdomain);
            this.LinearSystems = new ILinearSystem_v2[] { linearSystem };

            this.pcgAlgorithm = pcgAlgorithm;
            this.preconditionerFactory = preconditionerFactory;
            this.dofOrderer = dofOrderer;
        }

        public IReadOnlyList<ILinearSystem_v2> LinearSystems { get; }

        public IMatrix BuildGlobalMatrix(ISubdomain_v2 subdomain, IElementMatrixProvider elementMatrixProvider)
            => assembler.BuildGlobalMatrix(subdomain.DofOrdering, subdomain.Elements, elementMatrixProvider);

        public void Initialize()
        {
            // TODO: perhaps order dofs here
        }

        public void OrderDofs() { }
        //public void OrderDofs() => subdomain.DofOrdering = dofOrderer.OrderDofs(model, subdomain);

        /// <summary>
        /// Solves the linear system with PCG method. If the matrix has been modified, a new preconditioner will be computed.
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
            else linearSystem.Solution.Clear(); // In iterative algorithms we initialize the solution vector to 0.

            if (linearSystem.IsMatrixModified)
            {
                preconditioner = preconditionerFactory.CreatePreconditionerFor(linearSystem.Matrix);
                linearSystem.IsMatrixModified = false;
            }

            CGStatistics stats = pcgAlgorithm.Solve(linearSystem.Matrix, preconditioner, linearSystem.RhsVector,
                linearSystem.Solution, true); //TODO: This way, we don't know that x0=0, which will result in an extra b-A*0
        }

        //TODO: Create a method in Subdomain (or its DofOrderer) that exposes whether the dofs have changed.
        private bool HasSubdomainDofsChanged() => subdomain.DofOrdering.NumFreeDofs == linearSystem.Solution.Length;

        public class Builder
        {
            private MaxIterationsProvider maxIterationsProvider = new MaxIterationsProvider(1.0);

            public Builder() { }

            public IDofOrderer DofOrderer { get; set; } = new SimpleDofOrderer();

            public int MaxIterations
            {
                set => maxIterationsProvider = new MaxIterationsProvider(value);
            }

            public double MaxIterationsOverMatrixOrder
            {
                set => maxIterationsProvider = new MaxIterationsProvider(value);
            }

            public IPreconditionerFactory PreconditionerFactory { get; set; } = new JacobiPreconditioner.Factory();

            public double ResidualTolerance { get; set; } = 1E-6;

            public PcgSolver BuildSolver(IStructuralModel_v2 model)
            {
                return new PcgSolver(model, new PreconditionedConjugateGradient(maxIterationsProvider, ResidualTolerance),
                    PreconditionerFactory, DofOrderer);
            }
        }

        private class CsrSystem : LinearSystem_v2<CsrMatrix, Vector>
        {
            internal CsrSystem(ISubdomain_v2 subdomain) : base(subdomain) { }
            public override Vector CreateZeroVector() => Vector.CreateZero(Subdomain.DofOrdering.NumFreeDofs);
            public override void GetRhsFromSubdomain() => RhsVector = Vector.CreateFromArray(Subdomain.Forces, false);
        }
    }
}
