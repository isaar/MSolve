using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.SuiteSparse;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

//TODO: There are slight differences (e.g. in crack path) between this version ("no rebuild" and the older one "rebuild". 
//      Find the cause and fix it. I tried:
//      1) Requesting the columns of all taboo row indices: No change
//      2) Checking if the assembler builds all the requested columns: It does.
//      3) Checked the added columns of "no rebuild" and "rebuild": Same count. I bet that checking them individually will prove 
//          that they are identical. That needs to be done. Same needs to be done for removed columns.
//      4) Rebuilding the columns of all dofs in the modified elements, instead of only the requested ones: No change
//      5) Rebuilding all the elements of the model: Strangely this is closer to "no rebuilding" than "rebuilding whole K",
//          which makes me suspect that its less about the stiffnesses and more about the underlying matrices. Needs further
//          investigation. First check if the PartialMatrixColumns are identical with the columns of the whole rebuilt DOK.
//      6) Initializing PartialMatrixColumns as identity columns, like the DOK in "rebuild" version does: Surprisingly, this
//          produces almost identical results as the "rebuild" version, which revealed that the DOK initialization to unity was
//          incorrect in the first place. After fixing it, the "no rebuild" and "rebuild" versions produce ALMOST identical 
//          results. Currently it is not worth the time needed to find what causes this discrepancy, so it will be put on hold.
namespace ISAAR.MSolve.XFEM.Solvers
{
    /// <summary>
    /// In each iteration the factorization is updated with the added/removed columns. Only these columns are assembled, neither 
    /// the whole stiffness matrix nor the whole rhs vector. 
    /// </summary>
    class ReanalysisSolver : ISolver, IDisposable
    {
        private readonly TrackingExteriorCrackLSM crack;
        private readonly IReadOnlyList<XNode2D> fullyEnrichedNodes; // TODO: model must be passed in the constructor a parameter.
        private readonly Model2D model;

        private Vector prescribedNodalDisplacements;
        private Vector rhs;
        private CholeskySuiteSparse factorizedKff;
        private int counter;

        /// <summary>
        /// All nodes will be enriched with both Heaviside and crack tip functions to create the initial dof numbering. After 
        /// that, the enrichments will be cleared and only reapplied as needed by the crack propagation.
        /// </summary>
        public ReanalysisSolver(Model2D model, TrackingExteriorCrackLSM crack)
        {
            this.model = model;
            this.crack = crack;
            this.fullyEnrichedNodes = model.Nodes;
            Logger = new SolverLogger();
        }

        /// <summary>
        /// Only the provided nodes will be enriched with both Heaviside and crack tip functions to create the initial dof. After 
        /// that, the enrichments will be cleared and only reapplied as needed by the crack propagation. WARNING: if the crack 
        /// necessitates the enrichment of other nodes, an exception will be thrown.
        /// numbering.
        /// </summary>
        /// <param name="fullyEnrichedNodes">The nodes to be enriched with Heaviside and crack tip functions. Make sure that only
        ///     these nodes need to be enriched, otherwise an exception will be thrown.</param>
        public ReanalysisSolver(Model2D model, IReadOnlyList<XNode2D> fullyEnrichedNodes, TrackingExteriorCrackLSM crack)
        {
            this.model = model;
            this.crack = crack;
            this.fullyEnrichedNodes = fullyEnrichedNodes;
            Logger = new SolverLogger();
        }

        public IDofOrderer DofOrderer { get; protected set; }
        public SolverLogger Logger { get; }
        public Vector Solution { get; protected set; }

        public void Dispose()
        {
            if (factorizedKff != null) factorizedKff.Dispose();
        }

        public void Initialize()
        {
            var watch = new Stopwatch();
            watch.Start();

            // Enrich all applicable nodes, without evaluating the enrichment functions
            foreach (XNode2D node in fullyEnrichedNodes)
            {
                node.EnrichmentItems.Add(crack.CrackBodyEnrichment, null);
                node.EnrichmentItems.Add(crack.CrackTipEnrichments, null);
            }

            //DofOrderer = DofOrdererSeparate.Create(model);
            DofOrderer = InterleavedDofOrderer.Create(model);

            // Reorder and count the reordering time separately
            ReorderPatternSuperset(); //TODO: this actually increases fill-in

            // Clear all the possible enrichments. Only the required ones will be used as the crack propagates.
            foreach (XNode2D node in fullyEnrichedNodes) node.EnrichmentItems.Clear();

            watch.Stop();
            Logger.InitializationTime = watch.ElapsedMilliseconds;
        }

        public void Solve()
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

            var orderingAlgorithm = new OrderingAMD();
            (int[] permutationOldToNew, ReorderingStatistics stats) = pattern.Reorder(orderingAlgorithm);
            DofOrderer.ReorderUnconstrainedDofs(permutationOldToNew, false);
        }

        private void SolveFirstTime()
        {
            counter = 0;
            var watch = new Stopwatch();
            watch.Start();

            // Build the whole stiffness matrix for the first and last time.
            var assembler = new ReanalysisWholeAssembler();
            (DOKSymmetricColMajor Kff, CSRMatrix Kfc) = assembler.BuildGlobalMatrix(model.Elements, DofOrderer);

            /// The extended linear system is:
            /// [Kcc Kcf; Kuc Kff] * [uc; uf] = [Fc; Ff]
            /// where c are the standard constrained dofs, s are the standard free dofs, e are the enriched dofs and 
            /// where c are the standard constrained dofs, s are the standard free dofs, e are the enriched dofs and 
            /// f = Union(s,e) are both the dofs with unknown left hand side vectors: uu = [us; ue].
            /// To solve the system (for the unknowns ul):
            /// i) Kff * uf = Ff - Kfc * uc = Feff
            /// ii) uf = Kff \ Feff 
            Vector Ff = model.CalculateFreeForces(DofOrderer);
            prescribedNodalDisplacements = model.CalculateConstrainedDisplacements(DofOrderer);
            rhs = Ff - Kfc.MultiplyRight(prescribedNodalDisplacements); //TODO: directly multiply DOK*vector

            // Factorize the whole stiffness matrix for the first and last time.
            // WARNING: DO NOT use using(){} here. We want the unmanaged resource to persist for the lifetime of this object.
            factorizedKff = Kff.BuildSymmetricCSCMatrix(true).FactorCholesky(SuiteSparseOrdering.Natural); 
            Solution = factorizedKff.SolveLinearSystem(rhs);

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);

            //CheckSolutionAndEnforce(1.0);
            //CheckSolutionAndPrint(0.0);
        }

        //TODO: Try and compare the performance of first delete then add, modifying dofs in the order of their gloabl index. 
        //I think that modifying a matrix left to right is faster
        private void SolveByUpdating()
        {
            ++counter;
            var watch = new Stopwatch();
            watch.Start();

            // Group the new dofs
            IReadOnlyList<EnrichedDof> heavisideDofs = crack.CrackBodyEnrichment.Dofs;
            IReadOnlyList<EnrichedDof> tipDofs = crack.CrackTipEnrichments.Dofs;
            var colsToAdd = new HashSet<int>();
            foreach (var node in crack.CrackTipNodesNew)
            {
                foreach (var tipDof in tipDofs) colsToAdd.Add(DofOrderer.GetEnrichedDofOf(node, tipDof));
            }
            foreach (var node in crack.CrackBodyNodesNew)
            {
                foreach (var heavisideDof in heavisideDofs) colsToAdd.Add(DofOrderer.GetEnrichedDofOf(node, heavisideDof));
            }
            var tabooRows = new HashSet<int>(colsToAdd);

            // Delete old tip enriched Dofs. Should this be done before or after addition?
            foreach (var node in crack.CrackTipNodesOld)
            {
                foreach (var tipDof in tipDofs)
                {
                    int colIdx = DofOrderer.GetEnrichedDofOf(node, tipDof);
                    tabooRows.Add(colIdx); // They must also be excluded when adding new columns
                    factorizedKff.DeleteRow(colIdx);
                }
            }

            //WARNING: the next must happen after deleting old tip dofs and before adding new dofs
            TreatOtherModifiedDofs(heavisideDofs, colsToAdd, tabooRows);
            //Console.WriteLine($"{colsToAdd.Count} columns will be added");

            // Build only the stiffness matrix columns that must be added (even those I just removed). 
            // Also update RHS if needed.
            var elementsToRebuild = FindElementsToRebuild();
            var assembler = new ReanalysisPartialAssembler();
            PartialMatrixColumns changedStiffnessColumns = assembler.BuildGlobalMatrixColumns(
                elementsToRebuild, DofOrderer, colsToAdd, rhs, prescribedNodalDisplacements);

            // Add a column for each new dof. That column only contains entries corresponding to already active dofs and the 
            // new one. Adding the whole column is incorrect When a new column is added, that dof is removed from the set of 
            // remaining ones.
            foreach (int col in colsToAdd)
            {
                tabooRows.Remove(col); //It must be called before passing the remaining cols as taboo rows.
                SparseVector newColVector = changedStiffnessColumns.GetColumnWithoutRows(col, tabooRows);
                factorizedKff.AddRow(col, newColVector);
            }

            // Solve the system
            Solution = factorizedKff.SolveLinearSystem(rhs);

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);

            //CheckSolutionAndEnforce(1.0);
            //CheckSolutionAndPrint(0.0);
        }

        private HashSet<XContinuumElement2D> FindElementsToRebuild()
        {
            var elementsToRebuild = new HashSet<XContinuumElement2D>(crack.ElementsModified);

            /// The nodal support of <see cref="TrackingExteriorCrackLSM.CrackBodyNodesModified"/> includes elements with the 
            /// same stiffness matrix. However these matrices must be recomputed in order to build the global stiffness columns 
            /// corresponding to these nodes. 
            foreach (var node in crack.CrackBodyNodesNearModified)
            {
                foreach (var element in crack.Mesh.FindElementsWithNode(node)) elementsToRebuild.Add(element);
            }

            return elementsToRebuild;
            //return new HashSet<XContinuumElement2D>(model.Elements); // For debugging
        }

        private void TreatOtherModifiedDofs(IReadOnlyList<EnrichedDof> heavisideDofs, HashSet<int> colsToAdd,
            HashSet<int> tabooRows)
        {
            // Delete unmodified Heaviside dofs of nodes in elements with modified stiffness. 
            // Their stiffness will be readded later.
            // Same for previously Heaviside dofs with modified body level set.
            //TODO: Not sure, their stiffness can change. Investigate it.
            //TODO: Check if their stiffness is really modified. Otherwise it is a waste of time.
            foreach (var node in crack.CrackBodyNodesNearModified.Union(crack.CrackBodyNodesModified))
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
            (DOKSymmetricColMajor Kuu, CSRMatrix Kuc) = assembler.BuildGlobalMatrix(model.Elements, DofOrderer);
            Vector rhsNew = model.CalculateFreeForces(DofOrderer) 
                - Kuc.MultiplyRight(model.CalculateConstrainedDisplacements(DofOrderer));
            CholeskySuiteSparse factorization = Kuu.BuildSymmetricCSCMatrix(true).FactorCholesky(SuiteSparseOrdering.Natural);
            Vector solutionExpected = factorization.SolveLinearSystem(rhsNew);
            double error = (Solution - solutionExpected).Norm2() / solutionExpected.Norm2();
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
            string directory = @"C:\Users\Serafeim\Desktop\GRACM\Reanalysis_debugging\";
            string matrixPath = directory + "reanalysis_expected_matrix_" + counter + ".txt";
            string rhsPath = directory + "reanalysis_expected_rhs_" + counter + ".txt";
            string removedRowsPath = directory + "reanalysis_removed_rows_" + counter + ".txt";
            string addedRowsPath = directory + "reanalysis_added_rows_" + counter + ".txt";

            Console.WriteLine();
            Console.WriteLine("------------- DEBUG: reanalysis solver/ -------------");
            var assembler = new ReanalysisWholeAssembler();
            (DOKSymmetricColMajor Kuu, CSRMatrix Kuc) = assembler.BuildGlobalMatrix(model.Elements, DofOrderer);
            Vector rhsNew = model.CalculateFreeForces(DofOrderer)
                - Kuc.MultiplyRight(model.CalculateConstrainedDisplacements(DofOrderer));
            Vector solutionExpected;
            using (CholeskySuiteSparse factorization = Kuu.BuildSymmetricCSCMatrix(true).FactorCholesky(SuiteSparseOrdering.Natural))
            {
                solutionExpected = factorization.SolveLinearSystem(rhs);
            }
            double error = (Solution - solutionExpected).Norm2() / solutionExpected.Norm2();
            Console.Write($"Normalized error = {error}");

            if (error < tolerance)
            {
                Console.Write($". It is under the tolerance = {tolerance}.");
                Console.Write("Printing the expected matrix, rhs vector and altered dofs");

                // Expected matrix and rhs
                CoordinateTextFileSymmetricWriter.NumericFormat = new GeneralNumericFormat();
                (new CoordinateTextFileSymmetricWriter(Kuu)).WriteToFile(matrixPath);
                FullVectorWriter.NumericFormat = new GeneralNumericFormat();
                (new FullVectorWriter(rhs, true)).WriteToFile(rhsPath);

                // Modified dofs
                var removedRows = new List<int>();
                var addedRows = new List<int>();
                IReadOnlyList<EnrichedDof> heavisideDofs = crack.CrackBodyEnrichment.Dofs;
                IReadOnlyList<EnrichedDof> tipDofs = crack.CrackTipEnrichments.Dofs;
                foreach (var node in crack.CrackTipNodesOld)
                {
                    foreach (var tipDof in tipDofs) removedRows.Add(DofOrderer.GetEnrichedDofOf(node, tipDof));
                }
                foreach (var node in crack.CrackTipNodesNew)
                {
                    foreach (var tipDof in tipDofs) addedRows.Add(DofOrderer.GetEnrichedDofOf(node, tipDof));
                }
                foreach (var node in crack.CrackBodyNodesNew)
                {
                    foreach (var heavisideDof in heavisideDofs) addedRows.Add(DofOrderer.GetEnrichedDofOf(node, heavisideDof));
                }
                WriteList(removedRows, removedRowsPath);
                WriteList(addedRows, addedRowsPath);
            }

            Console.WriteLine();
            Console.WriteLine("------------- /DEBUG: reanalysis solver -------------");
            Console.WriteLine();
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
