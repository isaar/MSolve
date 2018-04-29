using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class ReanalysisSolver: SolverBase, IDisposable
    {
        private readonly IExteriorCrack crack;
        private readonly IReadOnlyList<XNode2D> fullyEnrichedNodes; // TODO: model must be passed in the constructor a parameter.
        private CholeskySuiteSparse factorizedKuu;

        /// <summary>
        /// All nodes will be enriched with both Heaviside and crack tip functions to create the initial dof numbering. After 
        /// that, the enrichments will be cleared and only reapplied as needed by the crack propagation.
        /// </summary>
        public ReanalysisSolver(Model2D model, IExteriorCrack crack) : base(model)
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
        public ReanalysisSolver(Model2D model, IReadOnlyList<XNode2D> fullyEnrichedNodes, IExteriorCrack crack) : base(model)
        {
            this.crack = crack;
            this.fullyEnrichedNodes = fullyEnrichedNodes;
        }

        public void Dispose()
        {
            if (factorizedKuu != null) factorizedKuu.Dispose();
        }

        public override void Initialize()
        {
            var watch = new Stopwatch();
            watch.Start();

            // Enrich all applicable nodes, without evaluating the enrichment functions
            foreach (XNode2D node in fullyEnrichedNodes) 
            {
                node.EnrichmentItems.Add(crack.CrackBodyEnrichment, null);
                node.EnrichmentItems.Add(crack.CrackTipEnrichments, null);
            }

            //DOFEnumerator = DOFEnumeratorSeparate.Create(model);
            DOFEnumerator = DOFEnumeratorInterleaved.Create(model);

            // Reorder and count the reordering time separately
            //
            // TODO: do that here 
            //

            // Clear all the possible enrichments. Only the required ones will be used as the crack propagates.
            foreach (XNode2D node in fullyEnrichedNodes) node.EnrichmentItems.Clear();

            // Build
            var assembler = new GlobalReanalysisAssembler();
            (DOKSymmetricColMajor Kuu, CSRMatrix Kuc) = assembler.BuildGlobalMatrix(model, DOFEnumerator);

            watch.Stop();
            Logger.InitializationTime = watch.ElapsedMilliseconds;
        }

        public override void Solve()
        {
            //TODO: it would be more clear if the SolveFirstTime() was done in Initialize(). However Initialize() may be called
            //      before the initial configuration of the model is built correclty. Another approach is to use a specialized 
            //      variation of QuasiStaticAnalysis and have control over when Initialize() is called.
            if (factorizedKuu == null) SolveFirstTime();
            else SolveByUpdating();
        }

        private void SolveFirstTime()
        {
            var watch = new Stopwatch();
            watch.Start();

            // TODO: the matrices must not be built and factorized in each iteration
            var assembler = new GlobalReanalysisAssembler();
            (DOKSymmetricColMajor Kuu, CSRMatrix Kuc) = assembler.BuildGlobalMatrix(model, DOFEnumerator);
            Vector rhs = CalcEffectiveRhs(Kuc);

            factorizedKuu = Kuu.BuildSymmetricCSCMatrix(true).FactorCholesky(); // DO NOT use using(){} here!
            Solution = factorizedKuu.SolveLinearSystem(rhs);

            //CheckSolution(0.0);

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);
        }

        private void SolveByUpdating()
        {
            var watch = new Stopwatch();
            watch.Start();

            // TODO: the matrices and vectors must not be built in each iteration
            var assembler = new GlobalReanalysisAssembler();
            (DOKSymmetricColMajor Kuu, CSRMatrix Kuc) = assembler.BuildGlobalMatrix(model, DOFEnumerator);
            Vector rhs = CalcEffectiveRhs(Kuc);

            // Delete old tip enriched DOFs
            IReadOnlyList<EnrichedDOF> tipDOFs = crack.CrackTipEnrichments.DOFs;
            foreach (var node in crack.CrackTipNodesOld)
            {
                foreach (var tipDOF in tipDOFs)
                {
                    int colIdx = DOFEnumerator.GetEnrichedDofOf(node, tipDOF);
                    factorizedKuu.DeleteRow(colIdx);
                }
            }

            // Add new tip DOFs
            foreach (var node in crack.CrackTipNodesNew)
            {
                foreach (var tipDOF in tipDOFs)
                {
                    int colIdx = DOFEnumerator.GetEnrichedDofOf(node, tipDOF);
                    SparseVector newCol = Kuu.BuildColumn(colIdx);
                    factorizedKuu.AddRow(colIdx, newCol);
                }
            }

            // Add new Heaviside DOFs
            IReadOnlyList<EnrichedDOF> heavisideDOFs = crack.CrackBodyEnrichment.DOFs;
            foreach (var node in crack.CrackBodyNodesNew)
            {
                foreach (var heavisideDOF in heavisideDOFs)
                {
                    int colIdx = DOFEnumerator.GetEnrichedDofOf(node, heavisideDOF);
                    SparseVector newCol = Kuu.BuildColumn(colIdx);
                    factorizedKuu.AddRow(colIdx, newCol);
                }
            }

            Solution = factorizedKuu.SolveLinearSystem(rhs);

            //CheckSolution(0.26);

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);
        }

        /// <summary>
        /// Can also set the computed solution vector to the expected, if the normalized error is under a specified tolerance.
        /// This is used to avoid or evaluate the sensitivity of other code components.
        /// </summary>
        /// <param name="solution"></param>
        /// <param name="enforceTolerance"></param>
        private void CheckSolution(double enforceTolerance)
        {
            Console.WriteLine();
            Console.WriteLine("------------- DEBUG: reanalysis solver/ -------------");
            var assembler = new GlobalReanalysisAssembler();
            (DOKSymmetricColMajor Kuu, CSRMatrix Kuc) = assembler.BuildGlobalMatrix(model, DOFEnumerator);
            Vector rhs = CalcEffectiveRhs(Kuc);
            CholeskySuiteSparse factorization = Kuu.BuildSymmetricCSCMatrix(true).FactorCholesky();
            Vector solutionExpected = factorization.SolveLinearSystem(rhs);
            double error = (Solution - solutionExpected).Norm2() / solutionExpected.Norm2();
            Console.Write($"Normalized error = {error}");
            if (error < enforceTolerance)
            {
                Console.Write($". It is under the tolerance = {enforceTolerance}.");
                Console.Write(" Setting the expected vector as solution. Also setting the factorized matrix to the correct one");
                Solution = solutionExpected;
                factorizedKuu.Dispose();
                factorizedKuu = factorization;
            }
            else
            {
                factorization.Dispose();
            }
            Console.WriteLine();
            Console.WriteLine("------------- /DEBUG: reanalysis solver -------------");
            Console.WriteLine();
        }
    }
}
