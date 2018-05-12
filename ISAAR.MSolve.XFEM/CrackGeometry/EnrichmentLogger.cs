using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Output.VTK;

//TODO: right now lsm.UpdateGeometry() is called before checking collapse, while lsm.UpdateEnrichments() afterwards. Thus there
//      is less 1 enrichment output file, compared to level set output files. Also if the geometry is updated before the analysis
//      starts there will not be any enrichment output files for those updates, whcih will cause a mismatch in the numbering.
namespace ISAAR.MSolve.XFEM.CrackGeometry
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

            // Log the Heaviside enriched nodes and the signs of their crack body level sets.
            writer.InitializeFile($"{outputDirectory}\\heaviside_nodes_{iteration}", true);
            var heavisideNodes = new Dictionary<ICartesianPoint2D, double>();
            foreach (var node in lsm.CrackBodyNodesAll)
            {
                double sign = Math.Sign(lsm.LevelSetsBody[node]);
                heavisideNodes.Add(node, sign);
            }
            writer.WriteScalarField("Heaviside_nodes", heavisideNodes);
            writer.CloseCurrentFile();

            // Log the tip enriched nodes and the signs of their crack body level sets.
            writer.InitializeFile($"{outputDirectory}\\tip_nodes_{iteration}", true);
            var tipNodes = new Dictionary<ICartesianPoint2D, double>();
            foreach (var node in lsm.CrackTipNodesNew)
            {
                double sign = Math.Sign(lsm.LevelSetsTip[node]);
                tipNodes.Add(node, sign);
            }
            writer.WriteScalarField("Tip_nodes", tipNodes);
            writer.CloseCurrentFile();

            // Log the nodes that belong to elements intersected by the crack, but are not enriched with Heaviside 
            writer.InitializeFile($"{outputDirectory}\\heaviside_rejected_nodes_{iteration}", true);
            var rejectedNodes = new Dictionary<ICartesianPoint2D, double>();
            foreach (var node in lsm.CrackBodyNodesRejected)
            {
                double sign = Math.Sign(lsm.LevelSetsBody[node]);
                rejectedNodes.Add(node, sign);
            }
            writer.WriteScalarField("Heaviside_rejected_nodes", rejectedNodes);
            writer.CloseCurrentFile();

            ++iteration;
        }
    }
}
