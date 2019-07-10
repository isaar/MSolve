using System;
using System.Collections.Generic;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Output.VTK;

//TODO: right now lsm.UpdateGeometry() is called before checking collapse, while lsm.UpdateEnrichments() afterwards. Thus there
//      is less 1 enrichment output file, compared to level set output files. Also if the geometry is updated before the analysis
//      starts there will not be any enrichment output files for those updates, whcih will cause a mismatch in the numbering.
namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.Logging
{
    public class EnrichmentLogger
    {
        private readonly TrackingExteriorCrackLSM lsm;
        private readonly XModel model;
        private readonly string outputDirectory;
        private int iteration;

        public EnrichmentLogger(XModel model, TrackingExteriorCrackLSM lsm, string outputDirectory)
        {
            this.model = model;
            this.lsm = lsm;
            this.outputDirectory = outputDirectory;
            iteration = 0;
        }

        public void Log()
        {
            // Log the new tip enriched nodes and the signs of their crack body level sets.
            using (var writer = new VtkPointWriter($"{outputDirectory}\\tip_nodes_new_{iteration}.vtk"))
            {
                var tipNodesNew = new Dictionary<CartesianPoint, double>();
                foreach (var node in lsm.CrackTipNodesNew[lsm.CrackTipEnrichments])
                {
                    double sign = Math.Sign(lsm.LevelSetsTip[node]);
                    tipNodesNew.Add(node, sign);
                }
                writer.WriteScalarField("Tip_nodes_new", tipNodesNew);
            }

            // Log the old tip enriched nodes and the signs of their crack body level sets.
            using (var writer = new VtkPointWriter($"{outputDirectory}\\tip_nodes_old_{iteration}.vtk"))
            {
                var tipNodesOld = new Dictionary<CartesianPoint, double>();
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
                    tipNodesOld.Add(new CartesianPoint(-0.01, -0.01), 0);
                }
                writer.WriteScalarField("Tip_nodes_old", tipNodesOld);
            }

            // Log all Heaviside enriched nodes and the signs of their crack body level sets.
            using (var writer = new VtkPointWriter($"{outputDirectory}\\heaviside_nodes_all_{iteration}.vtk"))
            {
                var heavisideNodesAll = new Dictionary<CartesianPoint, double>();
                foreach (var node in lsm.CrackBodyNodesAll[lsm.CrackBodyEnrichment])
                {
                    double sign = Math.Sign(lsm.LevelSetsBody[node]);
                    heavisideNodesAll.Add(node, sign);
                }
                writer.WriteScalarField("Heaviside_nodes_all", heavisideNodesAll);
            }

            // Log only the new Heaviside enriched nodes and the signs of their crack body level sets.
            using (var writer = new VtkPointWriter($"{outputDirectory}\\heaviside_nodes_new_{iteration}.vtk"))
            {
                var heavisideNodesNew = new Dictionary<CartesianPoint, double>();
                foreach (var node in lsm.CrackBodyNodesNew[lsm.CrackBodyEnrichment])
                {
                    double sign = Math.Sign(lsm.LevelSetsBody[node]);
                    heavisideNodesNew.Add(node, sign);
                }
                writer.WriteScalarField("Heaviside_nodes_new", heavisideNodesNew);
            }

            // Log the nodes that belong to elements intersected by the crack, but are not enriched with Heaviside 
            using (var writer = new VtkPointWriter($"{outputDirectory}\\heaviside_rejected_nodes_{iteration}.vtk"))
            {
                var rejectedNodes = new Dictionary<CartesianPoint, double>();
                foreach (var node in lsm.CrackBodyNodesRejected[lsm.CrackBodyEnrichment])
                {
                    double sign = Math.Sign(lsm.LevelSetsBody[node]);
                    rejectedNodes.Add(node, sign);
                }
                if (rejectedNodes.Count == 0) //a phony node just to get the Paraview reader working. TODO: find a more elegant solution.
                {
                    rejectedNodes.Add(new CartesianPoint(-0.01, -0.01), 0);
                }
                writer.WriteScalarField("Heaviside_rejected_nodes", rejectedNodes);
            }
                

            // Log unmodified Heaviside nodes of elements with at least one modified node
            using (var writer = new VtkPointWriter($"{outputDirectory}\\near_modified_heaviside_nodes_{iteration}.vtk"))
            {
                var nearModifiedHeavisideNodes = new Dictionary<CartesianPoint, double>();
                foreach (var node in lsm.CrackBodyNodesNearModified[lsm.CrackBodyEnrichment])
                {
                    double sign = Math.Sign(lsm.LevelSetsBody[node]);
                    nearModifiedHeavisideNodes.Add(node, sign);
                }
                if (nearModifiedHeavisideNodes.Count == 0) // a phony node just to get the Paraview reader working. TODO: find a more elegant solution.
                {
                    nearModifiedHeavisideNodes.Add(new CartesianPoint(-0.01, -0.01), 0);
                }
                writer.WriteScalarField("near_modified_heaviside_nodes", nearModifiedHeavisideNodes);
            }

            ++iteration;
        }
    }
}
