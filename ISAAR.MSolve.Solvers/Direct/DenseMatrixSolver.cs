using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
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
        private readonly bool isMatrixPositiveDefinite; //TODO: actually there should be 3 states: posDef, symmIndef, unsymm

        private bool factorizeInPlace = true;
        private bool mustInvert = true;
        private Matrix inverse;

        private DenseMatrixSolver(IStructuralModel model, IDofOrderer dofOrderer, bool isMatrixPositiveDefinite) :
            base(model, dofOrderer, new DenseMatrixAssembler(), "DenseMatrixSolver")
        {
            this.isMatrixPositiveDefinite = isMatrixPositiveDefinite;
        }

        public override Dictionary<int, IMatrix> BuildGlobalMatrices(IElementMatrixProvider elementMatrixProvider)
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

            var watch = new Stopwatch();
            watch.Start();
            Matrix matrix = assembler.BuildGlobalMatrix(subdomain.FreeDofOrdering, subdomain.Elements, elementMatrixProvider);
            watch.Stop();
            Logger.LogTaskDuration("Matrix assembly", watch.ElapsedMilliseconds);
            return new Dictionary<int, IMatrix> { { subdomain.ID, matrix } };
        }

        public override void HandleMatrixWillBeSet()
        {
            mustInvert = true;
            inverse = null;
        }

        public override void Initialize() {}

        public override void PreventFromOverwrittingSystemMatrices() => factorizeInPlace = false;

        /// <summary>
        /// Solves the linear system with back-forward substitution. If the matrix has been modified, it will be refactorized.
        /// </summary>
        public override void Solve()
        {
            var watch = new Stopwatch();
            if (linearSystem.SolutionConcrete == null) linearSystem.SolutionConcrete = linearSystem.CreateZeroVectorConcrete();
            //else linearSystem.Solution.Clear(); // no need to waste computational time on this in a direct solver

            // Factorization
            if (mustInvert)
            {
                watch.Start();
                if (isMatrixPositiveDefinite) inverse = linearSystem.Matrix.FactorCholesky(factorizeInPlace).Invert(true);
                else inverse = linearSystem.Matrix.FactorLU(factorizeInPlace).Invert(true);
                watch.Stop();
                Logger.LogTaskDuration("Matrix factorization", watch.ElapsedMilliseconds);
                watch.Reset();
                mustInvert = false;
            }

            // Substitutions
            watch.Start();
            inverse.MultiplyIntoResult(linearSystem.RhsConcrete, linearSystem.SolutionConcrete);
            watch.Stop();
            Logger.LogTaskDuration("Back/forward substitutions", watch.ElapsedMilliseconds);
            Logger.IncrementAnalysisStep();
        }

        protected override Matrix InverseSystemMatrixTimesOtherMatrix(IMatrixView otherMatrix)
        {
            var watch = new Stopwatch();
            if (mustInvert)
            {
                watch.Start();
                if (isMatrixPositiveDefinite) inverse = linearSystem.Matrix.FactorCholesky(factorizeInPlace).Invert(true);
                else inverse = linearSystem.Matrix.FactorLU(factorizeInPlace).Invert(true);
                watch.Stop();
                Logger.LogTaskDuration("Matrix factorization", watch.ElapsedMilliseconds);
                watch.Reset();
                mustInvert = false;
            }
            watch.Start();
            Matrix result = inverse.MultiplyRight(otherMatrix);
            watch.Stop();
            Logger.LogTaskDuration("Back/forward substitutions", watch.ElapsedMilliseconds);
            Logger.IncrementAnalysisStep();
            return result;
        }

        public class Builder: ISolverBuilder
        {
            public Builder() { }

            public IDofOrderer DofOrderer { get; set; }
                = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public bool IsMatrixPositiveDefinite { get; set; } = true;

            ISolver ISolverBuilder.BuildSolver(IStructuralModel model) => BuildSolver(model);

            public DenseMatrixSolver BuildSolver(IStructuralModel model)
                => new DenseMatrixSolver(model, DofOrderer, IsMatrixPositiveDefinite);
        }
    }
}
