using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms.CG;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;

//TODO: perhaps the user should choose the PCG settings himself and pass it. In this case, this should be named IterativeSolver.
//TODO: the maxIterations of PCG should be able to use the order of the matrix as a default value.
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
        private readonly ISubdomain subdomain;
        private readonly LinearSystem_v2<CsrMatrix, Vector> linearSystem;
        private readonly PreconditionedConjugateGradient pcgAlgorithm;
        private readonly IPreconditionerBuilder preconditionerBuilder;
        private FreeDofOrderer_v2 dofOrderer; //TODO: this should probably be accessed from the subdomain
        private IPreconditioner preconditioner;

        public PcgSolver(IStructuralModel model, int maxIterations, double residualTolerance,
            IPreconditionerBuilder preconditionerBuilder)
        {
            if (model.ISubdomainsDictionary.Count != 1) throw new InvalidSolverException(
                $"{name} can be used if there is only 1 subdomain");
            this.subdomain = model.ISubdomainsDictionary.First().Value;
            this.linearSystem = new LinearSystem_v2<CsrMatrix, Vector>(subdomain.ID);
            this.LinearSystems = new Dictionary<int, ILinearSystem_v2>(1) { { linearSystem.ID, linearSystem } };

            //TODO: resolve this weird dependency. The (Newmark)Analyzer needs the initial solution (which may be != 0, if it 
            // comes from a previous analysis) before the Solver has performed the first system solution. However, to initialize
            // it we need the rhs vector which is created when the user calls Model.ConnectDataStructures(). The correct would be
            // to only access the number of free dofs, but that also would be available after Model.ConnectDataStructures().
            this.linearSystem.Solution = Vector.CreateZero(subdomain.Forces.Length);

            this.pcgAlgorithm = new PreconditionedConjugateGradient(maxIterations, residualTolerance);
            this.preconditionerBuilder = preconditionerBuilder;
        }

        public IReadOnlyDictionary<int, ILinearSystem_v2> LinearSystems { get; }

        public IMatrix BuildGlobalMatrix(ISubdomain subdomain, IElementMatrixProvider elementMatrixProvider)
        {
            if (dofOrderer == null) dofOrderer = FreeDofOrderer_v2.CreateWithElementMajorFreeDofOrder(
                subdomain.ΙElementsDictionary.Values, subdomain.Constraints);
            return assembler.BuildGlobalMatrix(dofOrderer, subdomain.ΙElementsDictionary.Values, elementMatrixProvider);
        }

        public void Initialize()
        {
            // TODO: perhaps order dofs here
        }

        /// <summary>
        /// Solves the linear system with PCG method. If the matrix has been modified, a new preconditioner will be computed.
        /// </summary>
        public void Solve()
        {
            if (linearSystem.IsMatrixModified)
            {
                preconditioner = preconditionerBuilder.BuildPreconditioner(linearSystem.Matrix);
                linearSystem.IsMatrixModified = false;
            }

            CGStatistics stats = pcgAlgorithm.Solve(linearSystem.Matrix, preconditioner, linearSystem.RhsVector,
                linearSystem.Solution, false); //TODO: This way, we don't know that x0=0, which will result in an extra b-A*0
        }
    }
}
