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
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
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
        private readonly ICrackDescription crack;
        private readonly Model2D model;
        private readonly Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>> possibleEnrichments;

        private Vector prescribedNodalDisplacements;
        private Vector rhs;
        private CholeskySuiteSparse factorizedKff;
        private int iteration;

        /// <summary>
        /// All nodes will be enriched with both Heaviside and crack tip functions to create the initial dof numbering. After 
        /// that, the enrichments will be cleared and only reapplied as needed by the crack propagation.
        /// </summary>
        public ReanalysisSolver(Model2D model, ICrackDescription crack)
        {
            this.model = model;
            this.crack = crack;

            // Possible enrichments
            possibleEnrichments = new Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>>();
            foreach (var enrichment in crack.Enrichments)
            {
                possibleEnrichments.Add(enrichment, model.Nodes);
            }

            Logger = new SolverLogger("ReanalysisSolver");
        }

        /// <summary>
        /// Only the provided nodes will be enriched with all Heaviside and crack tip functions to create the initial dof. After 
        /// that, the enrichments will be cleared and only reapplied as needed by the crack propagation. WARNING: if the crack 
        /// necessitates the enrichment of other nodes, an exception will be thrown.
        /// numbering.
        /// </summary>
        /// <param name="fullyEnrichedNodes">The nodes to be enriched with Heaviside and crack tip functions. Make sure that only
        ///     these nodes need to be enriched, otherwise an exception will be thrown.</param>
        public ReanalysisSolver(Model2D model, TrackingExteriorCrackLSM crack, IReadOnlyList<XNode2D> fullyEnrichedNodes)
        {
            this.model = model;
            this.crack = crack;
            
            // Possible enrichments
            possibleEnrichments = new Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>>();
            possibleEnrichments.Add(crack.CrackBodyEnrichment, fullyEnrichedNodes);
            possibleEnrichments.Add(crack.CrackTipEnrichments, fullyEnrichedNodes);

            Logger = new SolverLogger("ReanalysisSolver");
        }

        public ReanalysisSolver(Model2D model, ICrackDescription crack,
            Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>> possibleEnrichments)
        {
            this.model = model;
            this.crack = crack;
            this.possibleEnrichments = possibleEnrichments;
            Logger = new SolverLogger("ReanalysisSolver");
        }

        public IDofOrderer DofOrderer { get; protected set; }
        public SolverLogger Logger { get; }
        public IVector Solution { get; protected set; }

        public void Dispose()
        {
            if (factorizedKff != null) factorizedKff.Dispose();
        }

        public void Initialize()
        {
            iteration = 0;
            var watch = new Stopwatch();
            watch.Start();

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

            watch.Stop();
            Logger.LogDuration(iteration, "AMD ordering", watch.ElapsedMilliseconds);
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
            ++iteration;
            var watch = new Stopwatch();

            // Build the whole stiffness matrix for the first and last time.
            watch.Start();
            var assembler = new ReanalysisWholeAssembler();
            (DokSymmetric Kff, DokRowMajor Kfc) = assembler.BuildGlobalMatrix(model.Elements, DofOrderer);

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

            watch.Stop();
            Logger.LogDuration(iteration, "linear system assembly", watch.ElapsedMilliseconds);

            // Factorize the whole stiffness matrix for the first and last time.
            // WARNING: DO NOT use using(){} here. We want the unmanaged resource to persist for the lifetime of this object.
            watch.Restart();
            factorizedKff = CholeskySuiteSparse.Factorize(Kff.BuildSymmetricCscMatrix(true), false);
            watch.Stop();
            Logger.LogDuration(iteration, "Cholesky factorization", watch.ElapsedMilliseconds);

            // Solve using the factorization
            watch.Restart();
            Solution = factorizedKff.SolveLinearSystem(rhs);
            watch.Stop();
            Logger.LogDuration(iteration, "back & forward substitution", watch.ElapsedMilliseconds);

            Logger.LogDofs(iteration, DofOrderer.NumStandardDofs + DofOrderer.NumEnrichedDofs);
            //CheckSolutionAndEnforce(1.0);
            //CheckSolutionAndPrint(0.0);
        }

        //TODO: Try and compare the performance of first delete then add, modifying dofs in the order of their gloabl index. 
        //I think that modifying a matrix left to right is faster
        private void SolveByUpdating()
        {
            ++iteration;
            var watch = new Stopwatch();
            long assemblyTime = 0, modifyTime = 0;

            Console.WriteLine($"Iteration {iteration}: ");
            Console.WriteLine($"Num standard dofs = " + DofOrderer.NumStandardDofs);


            // Group the new dofs
            watch.Start();
            var colsToAdd = new HashSet<int>();

            int numTipDofsNew = 0;
            foreach (var tipNodes in crack.CrackTipNodesNew)
            {
                IReadOnlyList<EnrichedDof> tipDofs = tipNodes.Key.Dofs;
                foreach (var node in tipNodes.Value)
                {
                    foreach (var tipDof in tipDofs)
                    {
                        colsToAdd.Add(DofOrderer.GetEnrichedDofOf(node, tipDof));
                        ++numTipDofsNew;
                    }
                }
            }
            Console.WriteLine("Num tip dofs added = " + numTipDofsNew);

            int numBodyDofsNew = 0;
            foreach (var bodyNodes in crack.CrackBodyNodesNew)
            {
                IReadOnlyList<EnrichedDof> heavisideDofs = bodyNodes.Key.Dofs;
                foreach (var node in bodyNodes.Value)
                {
                    foreach (var heavisideDof in heavisideDofs)
                    {
                        colsToAdd.Add(DofOrderer.GetEnrichedDofOf(node, heavisideDof));
                        ++numBodyDofsNew;
                    }
                }
            }
            var tabooRows = new HashSet<int>(colsToAdd);
            Console.WriteLine("Num Heaviside dofs added = " + numBodyDofsNew);


            // Delete old tip enriched Dofs. Should this be done before or after addition?
            int numTipDofsOld = 0;
            foreach (var tipNodes in crack.CrackTipNodesOld)
            {
                IReadOnlyList<EnrichedDof> tipDofs = tipNodes.Key.Dofs;
                foreach (var node in tipNodes.Value)
                {
                    foreach (var tipDof in tipDofs)
                    {
                        int colIdx = DofOrderer.GetEnrichedDofOf(node, tipDof);
                        tabooRows.Add(colIdx); // They must also be excluded when adding new columns
                        watch.Stop();
                        assemblyTime += watch.ElapsedMilliseconds;

                        // Delete column from cholesky factorization
                        watch.Restart();
                        factorizedKff.DeleteRow(colIdx);
                        watch.Stop();
                        modifyTime += watch.ElapsedMilliseconds;

                        ++numTipDofsOld;
                    }
                }
            }
            Console.WriteLine("Num tip dofs removed = " + numTipDofsOld);

            //WARNING: the next must happen after deleting old tip dofs and before adding new dofs
            TreatOtherModifiedDofs(colsToAdd, tabooRows, ref assemblyTime, ref modifyTime);
            //Console.WriteLine($"{colsToAdd.Count} columns will be added");

            // Build only the stiffness matrix columns that must be added (even those I just removed). 
            // Also update RHS if needed.
            watch.Restart();
            var elementsToRebuild = FindElementsToRebuild();
            var assembler = new ReanalysisPartialAssembler();
            PartialMatrixColumns changedStiffnessColumns = assembler.BuildGlobalMatrixColumns(
                elementsToRebuild, DofOrderer, colsToAdd, rhs, prescribedNodalDisplacements);
            watch.Stop();
            assemblyTime += watch.ElapsedMilliseconds;

            // Add a column for each new dof. That column only contains entries corresponding to already active dofs and the 
            // new one. Adding the whole column is incorrect When a new column is added, that dof is removed from the set of 
            // remaining ones.
            foreach (int col in colsToAdd)
            {
                // Create the column to add
                watch.Restart();
                tabooRows.Remove(col); //It must be called before passing the remaining cols as taboo rows.
                SparseVector newColVector = changedStiffnessColumns.GetColumnWithoutRows(col, tabooRows);
                watch.Stop();
                assemblyTime += watch.ElapsedMilliseconds;

                // Add the column to the cholesky factorization
                watch.Restart();
                factorizedKff.AddRow(col, newColVector);
                watch.Stop();
                modifyTime += watch.ElapsedMilliseconds;
            }

            Logger.LogDuration(iteration, "linear system assembly", assemblyTime);
            Logger.LogDuration(iteration, "Cholesky modifications", modifyTime);

            // Solve the system
            watch.Restart();
            Solution = factorizedKff.SolveLinearSystem(rhs);
            watch.Stop();
            Logger.LogDuration(iteration, "back & forward substitution", modifyTime);

            Logger.LogDofs(iteration, DofOrderer.NumStandardDofs + DofOrderer.NumEnrichedDofs);
            //CheckSolutionAndEnforce(1.0);
            //CheckSolutionAndPrint(0.0);
        }

        private HashSet<XContinuumElement2D> FindElementsToRebuild()
        {
            var elementsToRebuild = new HashSet<XContinuumElement2D>(crack.ElementsModified);

            /// The nodal support of <see cref="TrackingExteriorCrackLSM.crackBodyNodesModified"/> includes elements with the 
            /// same stiffness matrix. However these matrices must be recomputed in order to build the global stiffness columns 
            /// corresponding to these nodes.
            foreach (var bodyNodes in crack.CrackBodyNodesNearModified)
            {
                foreach (var node in bodyNodes.Value)
                {
                    foreach (var element in crack.Mesh.FindElementsWithNode(node)) elementsToRebuild.Add(element);
                }
            }

            return elementsToRebuild;
            //return new HashSet<XContinuumElement2D>(model.Elements); // For debugging
        }

        private void TreatOtherModifiedDofs(HashSet<int> colsToAdd,
            HashSet<int> tabooRows, ref long assemblyTime, ref long modifyTime)
        {
            var watch = new Stopwatch();
            watch.Start();
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
            watch.Stop();
            assemblyTime += watch.ElapsedMilliseconds;

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
                        watch.Restart();
                        int colIdx = DofOrderer.GetEnrichedDofOf(node, heavisideDof);
                        colsToAdd.Add(colIdx);
                        tabooRows.Add(colIdx); // They must also be excluded when adding new columns
                        watch.Stop();
                        assemblyTime += watch.ElapsedMilliseconds;

                        watch.Restart();
                        factorizedKff.DeleteRow(colIdx);
                        watch.Stop();
                        modifyTime += watch.ElapsedMilliseconds;
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
            Vector rhsNew = model.CalculateFreeForces(DofOrderer) 
                - Kuc.MultiplyRight(model.CalculateConstrainedDisplacements(DofOrderer));
            var factorization = CholeskySuiteSparse.Factorize(Kuu.BuildSymmetricCscMatrix(true), false);
            Vector solutionExpected = factorization.SolveLinearSystem(rhsNew);
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
            //string matrixPath = directory + "reanalysis_expected_matrix_" + iteration + ".txt";
            //string rhsPath = directory + "reanalysis_expected_rhs_" + iteration + ".txt";
            //string removedRowsPath = directory + "reanalysis_removed_rows_" + iteration + ".txt";
            //string addedRowsPath = directory + "reanalysis_added_rows_" + iteration + ".txt";

            //Console.WriteLine();
            //Console.WriteLine("------------- DEBUG: reanalysis solver/ -------------");
            //var assembler = new ReanalysisWholeAssembler();
            //(DOKSymmetricColMajor Kuu, DOKRowMajor Kuc) = assembler.BuildGlobalMatrix(model.Elements, DofOrderer);
            //Vector rhsNew = model.CalculateFreeForces(DofOrderer)
            //    - Kuc.MultiplyRight(model.CalculateConstrainedDisplacements(DofOrderer));
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
            //    IReadOnlyList<EnrichedDof> heavisideDofs = crack.DofsHeaviside;
            //    IReadOnlyList<EnrichedDof> tipDofs = crack.DofsTip;
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
