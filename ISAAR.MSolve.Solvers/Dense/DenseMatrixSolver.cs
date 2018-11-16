using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;

namespace ISAAR.MSolve.Solvers.Dense
{
    public class DenseMatrixSolver: ISolver_v2
    {
        private const string name = "SkylineSolver"; // for error messages
        private readonly DenseMatrixAssembler assembler = new DenseMatrixAssembler();
        private readonly ISubdomain subdomain;
        private readonly IDofOrderer dofOrderer;
        private readonly DenseSystem linearSystem;
        private IDofOrdering dofOrdering;
        private CholeskyFull factorizedMatrix;

        public DenseMatrixSolver(IStructuralModel model, IDofOrderer dofOrderer)
        {
            if (model.ISubdomainsDictionary.Count != 1) throw new InvalidSolverException(
                $"{name} can be used if there is only 1 subdomain");
            this.subdomain = model.ISubdomainsDictionary.First().Value;
            this.linearSystem = new DenseSystem(subdomain);
            this.LinearSystems = new Dictionary<int, ILinearSystem_v2>(1) { { subdomain.ID, linearSystem } };

            this.dofOrderer = dofOrderer;
        }

        public IReadOnlyDictionary<int, ILinearSystem_v2> LinearSystems { get; }

        public IMatrix BuildGlobalMatrix(ISubdomain subdomain, IElementMatrixProvider elementMatrixProvider)
        {
            if (dofOrdering == null) dofOrdering = dofOrderer.OrderDofs(subdomain);

            // DEBUG
            var writer = new FullMatrixWriter();
            writer.NumericFormat = new ExponentialFormat() { NumDecimalDigits = 2 };

            var dofOrderingSimple = (new SimpleDofOrderer()).OrderDofs(subdomain);
            Matrix simpleOrderK = (Matrix)assembler.BuildGlobalMatrix(
                dofOrderingSimple, subdomain.ΙElementsDictionary.Values, elementMatrixProvider);

            Console.WriteLine();
            Console.WriteLine("Global matrix with simple ordering");
            writer.WriteToConsole(simpleOrderK);

            var dofOrderingNodeMajor = (new NodeMajorDofOrderer()).OrderDofs(subdomain);
            Matrix nodeMajorK = (Matrix)assembler.BuildGlobalMatrix(
                dofOrderingNodeMajor, subdomain.ΙElementsDictionary.Values, elementMatrixProvider);

            Console.WriteLine();
            Console.WriteLine("Global matrix with node major ordering");
            writer.WriteToConsole(nodeMajorK);

            var permutationNodeMajorToSimple = new int[dofOrderingNodeMajor.NumFreeDofs];
            foreach ((INode node, DOFType dofType, int nodeMajorIdx) in dofOrderingNodeMajor.FreeDofs)
            {
                permutationNodeMajorToSimple[nodeMajorIdx] = dofOrderingSimple.FreeDofs[node, dofType];
            }

            Matrix reorderedK = nodeMajorK.Reorder(permutationNodeMajorToSimple, true);

            Console.WriteLine();
            Console.WriteLine("Global matrix with node major ordering, reordered to simple");
            writer.WriteToConsole(reorderedK);

            //Console.WriteLine("Existing global dof enumeration:");
            //Utilities.PrintDofOrder(subdomain);

            //Console.WriteLine();
            //Console.WriteLine("Using dof orderer:");
            //Console.WriteLine(dofOrderer.FreeDofs.ToString());

            //Console.WriteLine();
            //Console.WriteLine("Global matrix");
            //SkylineMatrix matrix = assembler.BuildGlobalMatrix(dofOrderer, subdomain.ΙElementsDictionary.Values, elementMatrixProvider);
            //var writer = new FullMatrixWriter();
            //writer.NumericFormat = new ExponentialFormat() { NumDecimalDigits = 2 };
            //writer.WriteToConsole(matrix);

            // END DEBUG

            return assembler.BuildGlobalMatrix(dofOrdering, subdomain.ΙElementsDictionary.Values, elementMatrixProvider);
        }

        //public IMatrix BuildGlobalMatrix(ISubdomain subdomain, IElementMatrixProvider elementMatrixProvider)
        //{
        //    return assembler.BuildGlobalMatrix(subdomain, elementMatrixProvider);
        //}

        public void Initialize()
        {
            // TODO: perhaps order dofs here
        }

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
                factorizedMatrix = linearSystem.Matrix.FactorCholesky();
                linearSystem.IsMatrixModified = false; //TODO: this is bad, since someone else might see it as unchanged. Better use observers.
                linearSystem.IsMatrixFactorized = true;
            }

            factorizedMatrix.SolveLinearSystem(linearSystem.RhsVector, linearSystem.Solution);
        }

        //TODO: Create a method in Subdomain (or its DofOrderer) that exposes whether the dofs have changed.
        /// <summary>
        /// The number of dofs might have been changed since the previous Solution vector had been created.
        /// </summary>
        private bool HasSubdomainDofsChanged() => subdomain.TotalDOFs == linearSystem.Solution.Length;

        public class Builder
        {
            public Builder() { }

            public IDofOrderer DofOrderer { get; set; } = new SimpleDofOrderer();

            public DenseMatrixSolver BuildSolver(IStructuralModel model)
                => new DenseMatrixSolver(model, DofOrderer);
        }

        private class DenseSystem : LinearSystem_v2<Matrix, Vector>
        {
            private readonly ISubdomain subdomain;
            internal DenseSystem(ISubdomain subdomain) : base(subdomain.ID) => this.subdomain = subdomain;
            public override Vector CreateZeroVector() => Vector.CreateZero(subdomain.TotalDOFs);
            public override void GetRhsFromSubdomain() => RhsVector = Vector.CreateFromArray(subdomain.Forces, false);
        }
    }
}
