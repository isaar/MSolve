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
    class ReanalysisSolver: SolverBase
    {
        private readonly IExteriorCrack crack;
        private IReadOnlyList<XNode2D> fullyEnrichedNodes; // TODO: model must be passed in the constructor a parameter.

        /// <summary>
        /// All nodes will be enriched with both Heaviside and crack tip functions to create the initial dof numbering. After 
        /// that, the enrichments will be cleared and only reapplied as needed by the crack propagation.
        /// </summary>
        public ReanalysisSolver(Model2D model, IExteriorCrack crack) : base(model)
        {
            this.crack = crack;
            this.fullyEnrichedNodes = null;
        }

        /// <summary>
        /// Only the provided nodes will be enriched with both Heaviside and crack tip functions to create the initial dof. After 
        /// that, the enrichments will be cleared and only reapplied as needed by the crack propagation. WARNING: if the crack 
        /// necessitates the enrichment of other nodes, an exception will be thrown.
        /// numbering.
        /// </summary>
        /// <param name="fullyEnrichedNodes">The nodes to be enriched with Heaviside and crack tip functions. Make sure that only
        ///     these nodes need to be enriched, otherwise an exception will be thrown.</param>
        public ReanalysisSolver(Model2D model, IReadOnlyList<XNode2D> fullyEnrichedNodes) : base(model)
        {
            this.fullyEnrichedNodes = fullyEnrichedNodes;
        }

        public override void Initialize()
        {
            var watch = new Stopwatch();
            watch.Start();

            if (fullyEnrichedNodes == null) fullyEnrichedNodes = model.Nodes;

            // Enrich all applicable nodes, without evaluating the enrichment functions
            foreach (XNode2D node in fullyEnrichedNodes) 
            {
                node.EnrichmentItems.Add(crack.CrackBodyEnrichment, null);
                node.EnrichmentItems.Add(crack.CrackTipEnrichments, null);
            }

            //DOFEnumerator = DOFEnumeratorSeparate.Create(model);
            DOFEnumerator = DOFEnumeratorInterleaved.Create(model);
            //TODO: reorder and count the reordering time separately
            
            // Clear all the possible enrichments. Only the required ones will be used as the crack propagates.
            foreach (XNode2D node in fullyEnrichedNodes) node.EnrichmentItems.Clear();

            watch.Stop();
            Logger.InitializationTime = watch.ElapsedMilliseconds;
        }

        public override void Solve()
        {
            var watch = new Stopwatch();
            watch.Start();

            // TODO: the matrices must not be built and factorized in each iteration
            var assembler = new GlobalReanalysisAssembler();
            (DOKSymmetricColMajor Kuu, CSRMatrix Kuc) = assembler.BuildGlobalMatrix(model, DOFEnumerator);
            Vector rhs = CalcEffectiveRhs(Kuc);

            using (CholeskySuiteSparse factorization = Kuu.BuildSymmetricCSCMatrix(true).FactorCholesky())
            {
                Solution = factorization.SolveLinearSystem(rhs);
            }

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);
        }
    }
}
