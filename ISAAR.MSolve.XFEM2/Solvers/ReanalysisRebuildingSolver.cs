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
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

//TODO: compare this to the reanalysis solver that only rebuilds certain columns to figure out the effect of certain dofs:
//      Heaviside dofs newa new body elements, old tip dofs, new tip dofs, etc. Also try varying the tip enrichment area.
namespace ISAAR.MSolve.XFEM.Solvers
{
    /// <summary>
    /// In each iteration the factorization is updated with the added/removed columns, but these are calculated by assembling 
    /// the whole stiffness matrix.
    /// </summary>
    class ReanalysisRebuildingSolver : SolverBase, IDisposable
    {
        private readonly ICrackDescription crack;
        private readonly Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>> possibleEnrichments;
        private CholeskySuiteSparse factorizedKff;
        private int counter;

        /// <summary>
        /// All nodes will be enriched with both Heaviside and crack tip functions to create the initial dof numbering. After 
        /// that, the enrichments will be cleared and only reapplied as needed by the crack propagation.
        /// </summary>
        public ReanalysisRebuildingSolver(Model2D model, ICrackDescription crack): base(model)
        {
            this.crack = crack;

            // Possible enrichments
            possibleEnrichments = new Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>>();
            foreach (var enrichment in crack.Enrichments)
            {
                possibleEnrichments.Add(enrichment, model.Nodes);
            }

            //Logger = new SolverLogger("ReanalysisSolver");
        }

        /// <summary>
        /// Only the provided nodes will be enriched with both Heaviside and crack tip functions to create the initial dof. After 
        /// that, the enrichments will be cleared and only reapplied as needed by the crack propagation. WARNING: if the crack 
        /// necessitates the enrichment of other nodes, an exception will be thrown.
        /// numbering.
        /// </summary>
        /// <param name="fullyEnrichedNodes">The nodes to be enriched with Heaviside and crack tip functions. Make sure that only
        ///     these nodes need to be enriched, otherwise an exception will be thrown.</param>
        public ReanalysisRebuildingSolver(Model2D model, IReadOnlyList<XNode2D> fullyEnrichedNodes, TrackingExteriorCrackLSM crack) : 
            base(model)
        {
            this.crack = crack;

            // Possible enrichments
            possibleEnrichments = new Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>>();
            possibleEnrichments.Add(crack.CrackBodyEnrichment, fullyEnrichedNodes);
            possibleEnrichments.Add(crack.CrackTipEnrichments, fullyEnrichedNodes);
        }

        public ReanalysisRebuildingSolver(Model2D model, ICrackDescription crack,
            Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>> possibleEnrichments) : base(model)
        {
            this.crack = crack;
            this.possibleEnrichments = possibleEnrichments;
        }

        public void Dispose()
        {
            if (factorizedKff != null) factorizedKff.Dispose();
        }

        public override void Initialize()
        {
            // Enrich all applicable nodes, without evaluating the enrichment functions
            foreach (var enrichmentNodes in possibleEnrichments)
            {
                foreach (XNode2D node in enrichmentNodes.Value) node.EnrichmentItems.Add(enrichmentNodes.Key, null);
            }

            //DofOrderer = DofOrdererSeparate.Create(model);
            DofOrderer = InterleavedDofOrderer.Create(model);

            // Reorder and count the reordering time separately
            ReorderPatternSuperset(); //TODO: this actually increases fill-in

            // Clear all the possible enrichments. Only the required ones will be used as the crack propagates.
            foreach (var enrichmentNodes in possibleEnrichments)
            {
                foreach (XNode2D node in enrichmentNodes.Value) node.EnrichmentItems.Clear();
            }
        }

        public override void Solve()
        {
            //TODO: it would be more clear if the SolveFirstTime() was done in Initialize(). However Initialize() may be called
            //      before the initial configuration of the model is built correctly. Another approach is to use a specialized 
            //      variation of QuasiStaticAnalysis and have control over when Initialize() is called.
            if (factorizedKff == null) SolveFirstTime();
            else SolveByUpdating();
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

            var orderingAlgorithm = new OrderingAmdSuiteSparse();
            (int[] permutation, bool oldToNew) = orderingAlgorithm.FindPermutation(pattern);
            DofOrderer.ReorderUnconstrainedDofs(permutation, oldToNew);
        }

        private void SolveFirstTime()
        {
            counter = 0;
            // TODO: the matrices must not be built and factorized in each iteration
            var assembler = new ReanalysisWholeAssembler();
            (DokSymmetric Kuu, DokRowMajor Kuc) = assembler.BuildGlobalMatrix(model.Elements, DofOrderer);
            Vector rhs = CalcEffectiveRhs(Kuc);

            factorizedKff = CholeskySuiteSparse.Factorize(Kuu.BuildSymmetricCscMatrix(true), false);
            Solution = factorizedKff.SolveLinearSystem(rhs);

            //CheckSolutionAndEnforce(1.0);
            CheckSolutionAndPrint(0.0);
        }

        //TODO: Try and compare the performance of first delete then add, modifying dofs in the order of their gloabl index. 
        //I think that modifying a matrix left to right is faster
        private void SolveByUpdating()
        {

            ++counter;
            // TODO: the matrices and vectors must not be built in each iteration
            var assembler = new ReanalysisWholeAssembler();
            (DokSymmetric Kff, DokRowMajor Kfc) = assembler.BuildGlobalMatrix(model.Elements, DofOrderer);
            Vector rhs = CalcEffectiveRhs(Kfc);

            // Group the new dofs
            var colsToAdd = new HashSet<int>();
            foreach (var tipNodes in crack.CrackTipNodesNew)
            {
                IReadOnlyList<EnrichedDof> tipDofs = tipNodes.Key.Dofs;
                foreach (var node in tipNodes.Value)
                {
                    foreach (var tipDof in tipDofs) colsToAdd.Add(DofOrderer.GetEnrichedDofOf(node, tipDof));
                }
            }
            foreach (var bodyNodes in crack.CrackBodyNodesNew)
            {
                IReadOnlyList<EnrichedDof> heavisideDofs = bodyNodes.Key.Dofs;
                foreach (var node in bodyNodes.Value)
                {
                    foreach (var heavisideDof in heavisideDofs) colsToAdd.Add(DofOrderer.GetEnrichedDofOf(node, heavisideDof));
                }
            }
            var tabooRows = new HashSet<int>(colsToAdd);

            // Delete old tip enriched Dofs. Should this be done before or after addition?
            foreach (var tipNodes in crack.CrackTipNodesOld)
            {
                IReadOnlyList<EnrichedDof> tipDofs = tipNodes.Key.Dofs;
                foreach (var node in tipNodes.Value)
                {
                    foreach (var tipDof in tipDofs)
                    {
                        int colIdx = DofOrderer.GetEnrichedDofOf(node, tipDof);
                        tabooRows.Add(colIdx); // They must also be excluded when adding new columns

                        // Delete column from cholesky factorization
                        factorizedKff.DeleteRow(colIdx);
                    }
                }
            }

            //WARNING: the next must happen after deleting old tip dofs and before adding new dofs
            TreatOtherModifiedDofs(colsToAdd, tabooRows);
            //Console.WriteLine($"{colsToAdd.Count} columns will be added");

            // Add a column for each new dof. That column only contains entries corresponding to already active dofs and the 
            // new one. Adding the whole column is incorrect When a new column is added, that dof is removed from the set of 
            // remaining ones.
            foreach (int col in colsToAdd)
            {
                tabooRows.Remove(col); //It must be called before passing the remaining cols as taboo rows.
                SparseVector newColVector = Kff.GetColumnWithoutRows(col, tabooRows);
                factorizedKff.AddRow(col, newColVector);
            }

            // Solve the system
            Solution = factorizedKff.SolveLinearSystem(rhs);

            //CheckSolutionAndEnforce(1.0);
            CheckSolutionAndPrint(0.0);
        }

        private void TreatOtherModifiedDofs(HashSet<int> colsToAdd, 
            HashSet<int> tabooRows)
        {
            var enrichments = new List<CrackBodyEnrichment2D>();
            var nodesToDelete = new List<ISet<XNode2D>>();
            foreach (var bodyNodes in crack.CrackBodyNodesModified)
            {
                enrichments.Add(bodyNodes.Key);
                nodesToDelete.Add(bodyNodes.Value);
            }
            foreach (var bodyNodes in crack.CrackBodyNodesNearModified)
            {
                enrichments.Add(bodyNodes.Key);
                nodesToDelete.Add(bodyNodes.Value);
            }

            // Delete unmodified Heaviside dofs of nodes in elements with modified stiffness. 
            // Their stiffness will be readded later.
            // Same for previously Heaviside dofs with modified body level set.
            //TODO: Not sure, their stiffness can change. Investigate it.
            //TODO: Check if their stiffness is really modified. Otherwise it is a waste of time.
            for (int i = 0; i < nodesToDelete.Count; ++i)
            {
                IReadOnlyList<EnrichedDof> heavisideDofs = enrichments[i].Dofs;
                foreach (var node in nodesToDelete[i])
                {
                    //Console.WriteLine($"Iteration {counter} - Near modified node: {node}");
                    foreach (var heavisideDof in heavisideDofs)
                    {
                        int colIdx = DofOrderer.GetEnrichedDofOf(node, heavisideDof);
                        colsToAdd.Add(colIdx);
                        tabooRows.Add(colIdx); // They must also be excluded when adding new columns

                        factorizedKff.DeleteRow(colIdx);
                    }
                }
            }
        }

        #region debugging
        /// <summary>
        /// Can also set the computed solution vector to the expected, if the normalized error is under a specified tolerance.
        /// This is used to avoid or evaluate the sensitivity of other code components.
        /// </summary>
        /// <param name="solution"></param>
        /// <param name="tolerance"></param>
        private void CheckSolutionAndEnforce(double tolerance)
        {
            Console.WriteLine();
            Console.WriteLine("------------- DEBUG: reanalysis solver/ -------------");
            var assembler = new ReanalysisWholeAssembler();
            (DokSymmetric Kuu, DokRowMajor Kuc) = assembler.BuildGlobalMatrix(model.Elements, DofOrderer);
            Vector rhs = CalcEffectiveRhs(Kuc);
            var factorization = CholeskySuiteSparse.Factorize(Kuu.BuildSymmetricCscMatrix(true), false);
            Vector solutionExpected = factorization.SolveLinearSystem(rhs);
            double error = Solution.Subtract(solutionExpected).Norm2() / solutionExpected.Norm2();
            Console.Write($"Normalized error = {error}");
            if (error < tolerance)
            {
                Console.Write($". It is under the tolerance = {tolerance}.");
                Console.Write(" Setting the expected vector as solution. Also setting the factorized matrix to the correct one");
                Solution = solutionExpected;
                factorizedKff.Dispose();
                factorizedKff = factorization;
            }
            else
            {
                factorization.Dispose();
            }
            Console.WriteLine();
            Console.WriteLine("------------- /DEBUG: reanalysis solver -------------");
            Console.WriteLine();
        }

        private void CheckSolutionAndPrint(double tolerance)
        {
            //string directory = @"C:\Users\Serafeim\Desktop\GRACM\Reanalysis_debugging\";
            //string matrixPath = directory + "reanalysis_expected_matrix_" + counter + ".txt";
            //string rhsPath = directory + "reanalysis_expected_rhs_" + counter + ".txt";
            //string removedRowsPath = directory + "reanalysis_removed_rows_" + counter + ".txt";
            //string addedRowsPath = directory + "reanalysis_added_rows_" + counter + ".txt";

            //Console.WriteLine();
            //Console.WriteLine("------------- DEBUG: reanalysis solver/ -------------");
            //var assembler = new ReanalysisWholeAssembler();
            //(DOKSymmetricColMajor Kuu, DOKRowMajor Kuc) = assembler.BuildGlobalMatrix(model.Elements, DofOrderer);
            //Vector rhs = CalcEffectiveRhs(Kuc);
            //Vector solutionExpected;
            //using (CholeskySuiteSparse factorization = Kuu.BuildSymmetricCSCMatrix(true).FactorCholesky())
            //{
            //    solutionExpected = factorization.SolveLinearSystem(rhs);
            //}
            //double error = (Solution - solutionExpected).Norm2() / solutionExpected.Norm2();
            //Console.Write($"Normalized error = {error}");

            //if (error < tolerance)
            //{
            //    Console.Write($". It is under the tolerance = {tolerance}.");
            //    Console.Write("Printing the expected matrix, rhs vector and altered dofs");

            //    // Expected matrix and rhs
            //    CoordinateTextFileSymmetricWriter.NumericFormat = new GeneralNumericFormat();
            //    (new CoordinateTextFileSymmetricWriter(Kuu)).WriteToFile(matrixPath);
            //    FullVectorWriter.NumericFormat = new GeneralNumericFormat();
            //    (new FullVectorWriter(rhs, true)).WriteToFile(rhsPath);

            //    // Modified dofs
            //    var removedRows = new List<int>();
            //    var addedRows = new List<int>();
            //    IReadOnlyList<EnrichedDof> heavisideDofs = crack.CrackBodyEnrichment.Dofs;
            //    IReadOnlyList<EnrichedDof> tipDofs = crack.CrackTipEnrichments.Dofs;
            //    foreach (var node in crack.CrackTipNodesOld)
            //    {
            //        foreach (var tipDof in tipDofs) removedRows.Add(DofOrderer.GetEnrichedDofOf(node, tipDof));
            //    }
            //    foreach (var node in crack.CrackTipNodesNew)
            //    {
            //        foreach (var tipDof in tipDofs) addedRows.Add(DofOrderer.GetEnrichedDofOf(node, tipDof));
            //    }
            //    foreach (var node in crack.CrackBodyNodesNew)
            //    {
            //        foreach (var heavisideDof in heavisideDofs) addedRows.Add(DofOrderer.GetEnrichedDofOf(node, heavisideDof));
            //    }
            //    WriteList(removedRows, removedRowsPath);
            //    WriteList(addedRows, addedRowsPath);
            //}

            //Console.WriteLine();
            //Console.WriteLine("------------- /DEBUG: reanalysis solver -------------");
            //Console.WriteLine();
        }

        private static void WriteList(List<int> dofs, string path)
        {
            using (var writer = new StreamWriter(path))
            {
                writer.WriteLine(dofs.Count);
                foreach (var dof in dofs) writer.Write(dof + " ");
            }
        }
        #endregion
    }
}
