using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.SuiteSparse;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

//TODO: time the assembly, reorder, factorization, solution, other separately
namespace ISAAR.MSolve.XFEM.Solvers
{
    // Only useful to assess the benefit of AMD
    class NoReanalysisSolver: SolverBase
    {
        private readonly IExteriorCrack crack;
        private readonly IReadOnlyList<XNode2D> fullyEnrichedNodes; // TODO: model must be passed in the constructor a parameter.
        private string dofOrdererName;
        private int[] permutationOldToNew;
        private IDofOrderer unorderedDofs;

        /// <summary>
        /// All nodes will be enriched with both Heaviside and crack tip functions to create the initial dof numbering. After 
        /// that, the enrichments will be cleared and only reapplied as needed by the crack propagation.
        /// </summary>
        public NoReanalysisSolver(Model2D model, IExteriorCrack crack) : base(model)
        {
            this.crack = crack;
            this.fullyEnrichedNodes = model.Nodes;
        }

        /// <summary>
        /// Only the provided nodes will be enriched with both Heaviside and crack tip functions to create the initial dof. After 
        /// that, the enrichments will be cleared and only reapplied as needed by the crack propagation. WARNING: if the crack 
        /// necessitates the enrichment of other nodes, an exception will be thrown.
        /// numbering.
        /// </summary>
        /// <param name="fullyEnrichedNodes">The nodes to be enriched with Heaviside and crack tip functions. Make sure that only
        ///     these nodes need to be enriched, otherwise an exception will be thrown.</param>
        public NoReanalysisSolver(Model2D model, IReadOnlyList<XNode2D> fullyEnrichedNodes, IExteriorCrack crack) : base(model)
        {
            this.crack = crack;
            this.fullyEnrichedNodes = fullyEnrichedNodes;
        }

        public override void Initialize()
        {
            // Enrich all applicable nodes, without evaluating the enrichment functions
            foreach (XNode2D node in fullyEnrichedNodes) 
            {
                node.EnrichmentItems.Add(crack.CrackBodyEnrichment, null);
                node.EnrichmentItems.Add(crack.CrackTipEnrichments, null);
            }

            //unorderedDofs = DofOrdererSeparate.Create(model);
            //dofOrdererName = "separate";
            unorderedDofs = InterleavedDofOrderer.Create(model);
            dofOrdererName = "interleaved";
            //unorderedDofs.WriteToConsole();

            // Keep a copy of the unordered dofs
            DofOrderer = unorderedDofs.DeepCopy();

            // Reorder and count the reordering time
            ReorderPatternSuperset();
            //DofOrderer.WriteToConsole();

            // Clear all the possible enrichments. Only the required ones will be used as the crack propagates.
            foreach (XNode2D node in fullyEnrichedNodes) node.EnrichmentItems.Clear();
        }

        public override void Solve()
        {
            var assembler = new ReanalysisWholeAssembler();
            //This is not posdef. Debugger showed empty columns. Also there are no enriched dofs!!!
            (DokSymmetric Kuu, DokRowMajor Kuc) = assembler.BuildGlobalMatrix(model.Elements, DofOrderer); 
            Vector rhs = CalcEffectiveRhs(Kuc);

            int nnzOrderedFactor;
            using (CholeskySuiteSparse factorization = Kuu.BuildSymmetricCscMatrix(true).FactorCholesky(SuiteSparseOrdering.Natural))
            {
                Solution = factorization.SolveLinearSystem(rhs);
                nnzOrderedFactor = factorization.NumNonZeros;
            }

            #region test
            int nnzUnorderedFactor = CheckSolutionWithUnordered();
            Console.WriteLine($"Initial ordering = {dofOrdererName}" 
                + $", reordering = AMD: nnz in factorization = {nnzOrderedFactor}");
            Console.WriteLine($"Initial ordering = {dofOrdererName}" 
                + $", reordering = none: nnz in factorization = {nnzUnorderedFactor}");
            Console.WriteLine();
            #endregion
        }

        private void ReorderPatternSuperset()
        {
            int order = DofOrderer.NumStandardDofs + DofOrderer.NumEnrichedDofs;
            var pattern = SparsityPatternSymmetric.CreateEmpty(order);

            // Could build the sparsity pattern during Dof enumeration?
            foreach (var element in model.Elements)
            {
                //TODO: what is the most efficient way to gather both? Perhaps the DofOrderer should do this 
                var standardDofs = DofOrderer.GetStandardDofsOf(element);
                var enrichedDofs = DofOrderer.GetEnrichedDofsOf(element);
                List<int> allDofs;
                if (standardDofs.Count > enrichedDofs.Count)
                {
                    allDofs = standardDofs;
                    allDofs.AddRange(enrichedDofs);
                }
                else
                {
                    allDofs = enrichedDofs;
                    allDofs.AddRange(standardDofs);
                }
                pattern.ConnectIndices(allDofs, false);
            }

            var orderingAlgorithm = new OrderingAmd();
            (int[] permutation, ReorderingStatistics stats) = orderingAlgorithm.FindPermutation(pattern);
            permutationOldToNew = permutation;
            DofOrderer.ReorderUnconstrainedDofs(permutationOldToNew, false);

        }

        /// <summary>
        /// Checks if without reordering, the solution is the same. Also returns the number of non zero entries in the 
        /// unordered factorized matrix.
        /// </summary>
        /// <returns></returns>
        private int CheckSolutionWithUnordered()
        {
            var assembler = new ReanalysisWholeAssembler();
            (DokSymmetric Kuu, DokRowMajor Kuc) = assembler.BuildGlobalMatrix(model.Elements, unorderedDofs);
            Vector Fu = model.CalculateFreeForces(unorderedDofs);
            Vector uc = model.CalculateConstrainedDisplacements(unorderedDofs);
            Vector rhs = Fu - Kuc.MultiplyRight(uc);
            Vector unorderSolution;
            int nnzUnorderedFactor;
            using (CholeskySuiteSparse factorization = Kuu.BuildSymmetricCscMatrix(true).FactorCholesky(SuiteSparseOrdering.Natural)) 
            {
                unorderSolution = factorization.SolveLinearSystem(rhs);
                nnzUnorderedFactor = factorization.NumNonZeros;
            }
            Vector solutionExpected = unorderSolution.Reorder(permutationOldToNew, false);
            double error = (Solution - solutionExpected).Norm2() / solutionExpected.Norm2();
            Console.WriteLine($"Normalized error compared to natural ordering = {error}");
            return nnzUnorderedFactor;
        }
    }
}
