using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class CholeskyAMDSolver: ISolver
    {
        protected readonly Model2D model;
        private string enumeratorName;
        private int iteration;

        public CholeskyAMDSolver(Model2D model)
        {
            this.model = model;
            Logger = new SolverLogger("CholeskyAMDSolver");
        }

        public IDofOrderer DofOrderer { get; private set; }
        public SolverLogger Logger { get; }
        public IVector Solution { get; private set; }

        public void Initialize()
        {
            iteration = 0;
        }

        public void Solve()
        {
            ++iteration;
            var watch = new Stopwatch();

            // Linear system assembly (part 1)
            watch.Start();
            var unorderedDofs = InterleavedDofOrderer.Create(model);
            enumeratorName = "interleaved";
            //var unorderedDofs = DofOrdererSeparate.Create(model);
            //enumeratorName = "separate";
            watch.Stop();
            long assemblyTime = watch.ElapsedMilliseconds;

            // Reordering
            watch.Restart();
            DofOrderer = unorderedDofs.DeepCopy(); //TODO: this just wastes time, doesn't it? I think it was for testing purposes.
            Reorder();
            watch.Stop();
            Logger.LogDuration(iteration, "AMD ordering", watch.ElapsedMilliseconds);

            // Linear system assembly (part 2)
            watch.Restart();
            var assembler = new GlobalDOKAssembler();
            (DokSymmetric Kuu, DokRowMajor Kuc) = assembler.BuildGlobalMatrix(model, DofOrderer);
            Vector rhs = CalcEffectiveRhs(Kuc);
            watch.Stop();
            assemblyTime += watch.ElapsedMilliseconds;
            Logger.LogDuration(iteration, "linear system assembly", assemblyTime);

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

            // Linear system solution
            watch.Restart();
            using (var factor = CholeskySuiteSparse.Factorize(Kuu.BuildSymmetricCscMatrix(true), true))
            {
                watch.Stop();
                Logger.LogDuration(iteration, "Cholesky factorization", watch.ElapsedMilliseconds);

                watch.Restart();
                Solution = factor.SolveLinearSystem(rhs);
                watch.Stop();
                Logger.LogDuration(iteration, "back & forward substitution", watch.ElapsedMilliseconds);
                //Console.WriteLine($"Ordering {enumeratorName} + AMD, nnz after factorization = {factorization.NumNonZeros}");
            }
            //SolveWithoutReordering(unorderedDofs);

            Logger.LogDofs(iteration, DofOrderer.NumStandardDofs + DofOrderer.NumEnrichedDofs);
        }

        private void Reorder()
        {
            var orderingAlgorithm = new OrderingAmdSuiteSparse();

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
            (int[] permutation, bool oldToNew) = orderingAlgorithm.FindPermutation(pattern);

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

            DofOrderer.ReorderUnconstrainedDofs(permutation, oldToNew);
        }

        private Vector CalcEffectiveRhs(DokRowMajor globalUnconstrainedConstrained)
        {
            Vector Fu = model.CalculateFreeForces(DofOrderer);
            Vector uc = model.CalculateConstrainedDisplacements(DofOrderer);
            Vector Feff = Fu - globalUnconstrainedConstrained.MultiplyRight(uc);
            return Feff;
        }

        private void SolveWithoutReordering(IDofOrderer unordered)
        {
            var assembler = new GlobalDOKAssembler();
            (DokSymmetric Kuu, DokRowMajor Kuc) = assembler.BuildGlobalMatrix(model, unordered);
            using (var factor = CholeskySuiteSparse.Factorize(Kuu.BuildSymmetricCscMatrix(true), true))
            {
                //Solution = factorization.SolveLinearSystem(rhs);
                Console.WriteLine($"Ordering {enumeratorName} unordered, nnz after factorization = {factor.NumNonZerosUpper}");
            }
        }
    }
}
