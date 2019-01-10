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
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.Direct
{
    public class DenseMatrixSolver: ISolver_v2
    {
        private const string name = "SkylineSolver"; // for error messages
        private readonly DenseMatrixAssembler assembler = new DenseMatrixAssembler();
        private readonly IDofOrderer dofOrderer;
        private readonly IStructuralModel_v2 model;
        private readonly ISubdomain_v2 subdomain;
        private readonly DenseSystem linearSystem;

        private bool mustFactorize = true;
        private CholeskyFull factorizedMatrix;

        public DenseMatrixSolver(IStructuralModel_v2 model, IDofOrderer dofOrderer)
        {
            if (model.Subdomains.Count != 1) throw new InvalidSolverException(
                $"{name} can be used if there is only 1 subdomain");
            this.model = model;
            subdomain = model.Subdomains[0];

            linearSystem = new DenseSystem(subdomain);
            LinearSystems = new Dictionary<int, ILinearSystem_v2>() { { subdomain.ID, linearSystem } };
            linearSystem.MatrixObservers.Add(this);

            this.dofOrderer = dofOrderer;
        }

        public IReadOnlyDictionary<int, ILinearSystem_v2> LinearSystems { get; }

        public IMatrix BuildGlobalMatrix(ISubdomain_v2 subdomain, IElementMatrixProvider_v2 elementMatrixProvider)
        {
            #region Code to facilitate debugging
            //var writer = new FullMatrixWriter();
            //writer.NumericFormat = new ExponentialFormat() { NumDecimalDigits = 2 };

            //var dofOrderingSimple = (new SimpleDofOrderer()).OrderDofs(model, subdomain);
            //Matrix simpleOrderK = (Matrix)assembler.BuildGlobalMatrix(
            //    dofOrderingSimple, subdomain.Elements, elementMatrixProvider);

            //Console.WriteLine();
            //Console.WriteLine("Global matrix with simple ordering");
            //writer.WriteToConsole(simpleOrderK);

            //var dofOrderingNodeMajor = (new NodeMajorDofOrderer()).OrderDofs(model, subdomain);
            //Matrix nodeMajorK = (Matrix)assembler.BuildGlobalMatrix(
            //    dofOrderingNodeMajor, subdomain.Elements, elementMatrixProvider);

            //Console.WriteLine();
            //Console.WriteLine("Global matrix with node major ordering");
            //writer.WriteToConsole(nodeMajorK);

            //var permutationNodeMajorToSimple = new int[dofOrderingNodeMajor.NumFreeDofs];
            //foreach ((INode node, DOFType dofType, int nodeMajorIdx) in dofOrderingNodeMajor.FreeDofs)
            //{
            //    permutationNodeMajorToSimple[nodeMajorIdx] = dofOrderingSimple.FreeDofs[node, dofType];
            //}

            //Matrix reorderedK = nodeMajorK.Reorder(permutationNodeMajorToSimple, true);

            //Console.WriteLine();
            //Console.WriteLine("Global matrix with node major ordering, reordered to simple");
            //writer.WriteToConsole(reorderedK);

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
            #endregion

            return assembler.BuildGlobalMatrix(subdomain.DofOrdering, subdomain.Elements, 
                elementMatrixProvider);
        }

        public void Initialize() {}

        public void OnMatrixSetting()
        {
            mustFactorize = true;
            factorizedMatrix = null;
        }

        public void OrderDofsAndClearLinearSystem()
        {
            assembler.OnDofOrderingModified();
            model.GlobalDofOrdering = dofOrderer.OrderDofs(model);
            OnMatrixSetting();
            linearSystem.Clear();
            linearSystem.Size = model.GlobalDofOrdering.SubdomainDofOrderings[subdomain].NumFreeDofs;
        }

        /// <summary>
        /// Solves the linear system with back-forward substitution. If the matrix has been modified, it will be refactorized.
        /// </summary>
        public void Solve()
        {

            if (linearSystem.Solution == null) linearSystem.Solution = linearSystem.CreateZeroVector();
            //else linearSystem.Solution.Clear(); // no need to waste computational time on this in a direct solver

            if (mustFactorize)
            {
                factorizedMatrix = linearSystem.Matrix.FactorCholesky();
                mustFactorize = false;
                linearSystem.IsMatrixOverwrittenBySolver = true;
            }

            factorizedMatrix.SolveLinearSystem(linearSystem.RhsVector, linearSystem.Solution);
        }

        public class Builder
        {
            public Builder() { }

            public IDofOrderer DofOrderer { get; set; }
                = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public DenseMatrixSolver BuildSolver(IStructuralModel_v2 model)
                => new DenseMatrixSolver(model, DofOrderer);
        }

        private class DenseSystem : LinearSystem_v2<Matrix, Vector>
        {
            internal DenseSystem(ISubdomain_v2 subdomain) : base(subdomain) { }
            public override void GetRhsFromSubdomain() => RhsVector = Subdomain.Forces;
            internal override Vector CreateZeroVector() => Vector.CreateZero(Size);
        }
    }
}
