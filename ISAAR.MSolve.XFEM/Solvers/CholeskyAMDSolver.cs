using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.SuiteSparse;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class CholeskyAMDSolver: SolverBase
    {
        private string enumeratorName;

        public CholeskyAMDSolver(Model2D model) : base(model)
        { }

        public override void Solve()
        {
            var watch = new Stopwatch();
            watch.Start();

            var unorderedDofs = InterleavedDofOrderer.Create(model);
            enumeratorName = "interleaved";
            //var unorderedDofs = DofOrdererSeparate.Create(model);
            //enumeratorName = "separate";

            DofOrderer = unorderedDofs.DeepCopy();
            Reorder();

            var assembler = new GlobalDOKAssembler();
            (DOKSymmetricColMajor Kuu, DOKRowMajor Kuc) = assembler.BuildGlobalMatrix(model, DofOrderer);
            Vector rhs = CalcEffectiveRhs(Kuc);

            #region debug
            //// Matrix
            //(double[] values, int[] rowIndices, int[] colOffsets) = Kuu.BuildSymmetricCSCArrays(true);
            //var csc = CSCMatrix.CreateFromArrays(Kuu.NumRows, Kuu.NumColumns, values, rowIndices, colOffsets, false);
            //(new RawArraysWriter(csc)).WriteToMultipleFiles(matrixPath, true);

            //// Rhs
            //FullVectorWriter.NumericFormat = new GeneralNumericFormat();
            //(new FullVectorWriter(rhs, true)).WriteToFile(rhsPath);

            //// Solution
            //Vector expectedSolution = SolveWithSkyline();
            //(new FullVectorWriter(expectedSolution, true)).WriteToFile(solutionPath);
            #endregion

            using (CholeskySuiteSparse factorization = Kuu.BuildSymmetricCSCMatrix(true).FactorCholesky(SuiteSparseOrdering.Natural))
            {
                Solution = factorization.SolveLinearSystem(rhs);
                Console.WriteLine($"Ordering {enumeratorName} + AMD, nnz after factorization = {factorization.NumNonZeros}");
            }
            SolveWithoutReordering(unorderedDofs);

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);
        }

        private void Reorder()
        {
            var orderingAlgorithm = new OrderingAMD();

            int order = DofOrderer.NumStandardDofs + DofOrderer.NumEnrichedDofs;
            var pattern = SparsityPatternSymmetricColMajor.CreateEmpty(order);
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
            (int[] permutationOldToNew, ReorderingStatistics stats) = pattern.Reorder(orderingAlgorithm);

            #region DEBUG
            //var assembler = new GlobalDOKAssembler();
            //(DOKSymmetricColMajor Kuu, DOKRowMajor Kuc) = assembler.BuildGlobalMatrix(model, DofOrderer);
            //(int[] permutation2, ReorderingStatistics stats2) = Kuu.Reorder(orderingAlgorithm);
            //var comparer = new Comparer();
            //bool sameReordering = comparer.AreEqual(permutationOldToNew, permutation2);
            //if (sameReordering) Console.WriteLine("Sparsity pattern and DOK produce the same reordering");
            //else Console.WriteLine("Sparsity pattern and DOK produce differenr reorderings");
            //Console.WriteLine($"Ordering {enumeratorName} + AMD, factor nnz predicted = {stats.FactorizedNumNonZeros}");
            #endregion

            DofOrderer.ReorderUnconstrainedDofs(permutationOldToNew, false);
        }

        private void SolveWithoutReordering(IDofOrderer unordered)
        {
            var assembler = new GlobalDOKAssembler();
            (DOKSymmetricColMajor Kuu, DOKRowMajor Kuc) = assembler.BuildGlobalMatrix(model, unordered);
            using (CholeskySuiteSparse factorization = Kuu.BuildSymmetricCSCMatrix(true).FactorCholesky(SuiteSparseOrdering.Natural))
            {
                //Solution = factorization.SolveLinearSystem(rhs);
                Console.WriteLine($"Ordering {enumeratorName} unordered, nnz after factorization = {factorization.NumNonZeros}");
            }
        }


    }
}
