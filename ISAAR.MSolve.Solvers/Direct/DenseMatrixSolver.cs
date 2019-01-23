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
    /// Direct solver for models with only 1 subdomain. Uses Cholesky factorization on symmetric positive definite matrices
    /// stored in full format. Its purpose is mainly for testing, since it is inefficient for large linear systems resulting 
    /// from FEM .
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class DenseMatrixSolver: SingleSubdomainSolverBase<Matrix>
    {
        private bool factorizeInPlace = true;
        private bool mustFactorize = true;
        private CholeskyFull factorizedMatrix;

        private DenseMatrixSolver(IStructuralModel_v2 model, IDofOrderer dofOrderer):
            base(model, dofOrderer, new DenseMatrixAssembler(), "DenseMatrixSolver")
        { }

        public override IMatrix BuildGlobalMatrix(ISubdomain_v2 subdomain, IElementMatrixProvider_v2 elementMatrixProvider)
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

        public override void Initialize() {}

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
                factorizedMatrix = linearSystem.Matrix.FactorCholesky(factorizeInPlace);
                mustFactorize = false;
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
    }
}
