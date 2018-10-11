using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Output.VTK;

//TODO: right now lsm.UpdateGeometry() is called before checking collapse, while lsm.UpdateEnrichments() afterwards. Thus there
//      is less 1 enrichment output file, compared to level set output files. Also if the geometry is updated before the analysis
//      starts there will not be any enrichment output files for those updates, whcih will cause a mismatch in the numbering.
namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.Logging
{
    class EnrichmentLogger
    {
        private readonly TrackingExteriorCrackLSM lsm;
        private readonly Model2D model;
        private readonly string outputDirectory;
        private int iteration;

        public EnrichmentLogger(Model2D model, TrackingExteriorCrackLSM lsm, string outputDirectory)
        {
            this.model = model;
            this.lsm = lsm;
            this.outputDirectory = outputDirectory;
            iteration = 0;
        }

        public void Log()
        {
            var writer = new VTKPointWriter();
            
            // Log the new tip enriched nodes and the signs of their crack body level sets.
            writer.InitializeFile($"{outputDirectory}\\tip_nodes_new_{iteration}", true);
            var tipNodesNew = new Dictionary<ICartesianPoint2D, double>();
            foreach (var node in lsm.CrackTipNodesNew[lsm.CrackTipEnrichments])
            {
                double sign = Math.Sign(lsm.LevelSetsTip[node]);
                tipNodesNew.Add(node, sign);
            }
            writer.WriteScalarField("Tip_nodes_new", tipNodesNew);
            writer.CloseCurrentFile();

            // Log the old tip enriched nodes and the signs of their crack body level sets.
            writer.InitializeFile($"{outputDirectory}\\tip_nodes_old_{iteration}", true);
            var tipNodesOld = new Dictionary<ICartesianPoint2D, double>();
            if (iteration > 0) 
            {
                foreach (var node in lsm.CrackTipNodesOld[lsm.CrackTipEnrichments])
                {
                    double sign = Math.Sign(lsm.LevelSetsTip[node]);
                    tipNodesOld.Add(node, sign);
                }
            }
            else // else a phony node just to get the Paraview reader working. TODO: find a more elegant solution.
            {
                tipNodesOld.Add(new CartesianPoint2D(-0.01, -0.01), 0);
            }
            writer.WriteScalarField("Tip_nodes_old", tipNodesOld);
            writer.CloseCurrentFile();

            // Log all Heaviside enriched nodes and the signs of their crack body level sets.
            writer.InitializeFile($"{outputDirectory}\\heaviside_nodes_all_{iteration}", true);
            var heavisideNodesAll = new Dictionary<ICartesianPoint2D, double>();
            foreach (var node in lsm.CrackBodyNodesAll[lsm.CrackBodyEnrichment])
            {
                double sign = Math.Sign(lsm.LevelSetsBody[node]);
                heavisideNodesAll.Add(node, sign);
            }
            writer.WriteScalarField("Heaviside_nodes_all", heavisideNodesAll);
            writer.CloseCurrentFile();

            // Log only the new Heaviside enriched nodes and the signs of their crack body level sets.
            writer.InitializeFile($"{outputDirectory}\\heaviside_nodes_new_{iteration}", true);
            var heavisideNodesNew = new Dictionary<ICartesianPoint2D, double>();
            foreach (var node in lsm.CrackBodyNodesNew[lsm.CrackBodyEnrichment])
            {
                double sign = Math.Sign(lsm.LevelSetsBody[node]);
                heavisideNodesNew.Add(node, sign);
            }
            writer.WriteScalarField("Heaviside_nodes_new", heavisideNodesNew);
            writer.CloseCurrentFile();

            // Log the nodes that belong to elements intersected by the crack, but are not enriched with Heaviside 
            writer.InitializeFile($"{outputDirectory}\\heaviside_rejected_nodes_{iteration}", true);
            var rejectedNodes = new Dictionary<ICartesianPoint2D, double>();
            foreach (var node in lsm.CrackBodyNodesRejected[lsm.CrackBodyEnrichment])
            {
                double sign = Math.Sign(lsm.LevelSetsBody[node]);
                rejectedNodes.Add(node, sign);
            }
            if (rejectedNodes.Count == 0) //a phony node just to get the Paraview reader working. TODO: find a more elegant solution.
            {
                rejectedNodes.Add(new CartesianPoint2D(-0.01, -0.01), 0);
            }
            writer.WriteScalarField("Heaviside_rejected_nodes", rejectedNodes);
            writer.CloseCurrentFile();


            // Log unmodified Heaviside nodes of elements with at least one modified node
            writer.InitializeFile($"{outputDirectory}\\near_modified_heaviside_nodes_{iteration}", true);
            var nearModifiedHeavisideNodes = new Dictionary<ICartesianPoint2D, double>();
            foreach (var node in lsm.CrackBodyNodesNearModified[lsm.CrackBodyEnrichment])
            {
                double sign = Math.Sign(lsm.LevelSetsBody[node]);
                nearModifiedHeavisideNodes.Add(node, sign);
            }
            if (nearModifiedHeavisideNodes.Count == 0) // a phony node just to get the Paraview reader working. TODO: find a more elegant solution.
            {
                nearModifiedHeavisideNodes.Add(new CartesianPoint2D(-0.01, -0.01), 0);
            }
            writer.WriteScalarField("near_modified_heaviside_nodes", nearModifiedHeavisideNodes);
            writer.CloseCurrentFile();

            ++iteration;

        }
    }
}
