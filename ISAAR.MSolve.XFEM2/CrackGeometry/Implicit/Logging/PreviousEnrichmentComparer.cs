using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.Logging
{
    class PreviousLevelSetComparer
    {
        private readonly TrackingExteriorCrackLSM lsm;
        private readonly string outputDirectory;
        private int iteration;
        private Dictionary<XNode2D, double> previousBodyLevelSets;
        private Dictionary<XNode2D, double> previousTipLevelSets;

        public PreviousLevelSetComparer(TrackingExteriorCrackLSM lsm, string outputDirectory)
        {
            this.lsm = lsm;
            this.outputDirectory = outputDirectory;
            iteration = 0;
        }

        public void Log()
        {
            if (iteration > 0)
            {
                using (var writer = new StreamWriter($"{outputDirectory}\\level_set_comparison_{iteration}.txt", false))
                {
                    writer.WriteLine("Old Heaviside nodes, body level set:");
                    foreach (var node in lsm.CrackBodyNodesAll[lsm.CrackBodyEnrichment].
                        Except(lsm.CrackBodyNodesNew[lsm.CrackBodyEnrichment]))
                    {
                        double previous = previousBodyLevelSets[node];
                        double current = lsm.LevelSetsBody[node];
                        string comparison = (previous == current) ? "SAME" : "DIFFERENT";
                        writer.WriteLine(
                            $"{node}: {comparison} - previous body level set = {previous} - current body level set = {current}");
                    }
                    writer.WriteLine();

                    writer.WriteLine("New Heaviside nodes, body level set:");
                    foreach (var node in lsm.CrackBodyNodesNew[lsm.CrackBodyEnrichment])
                    {
                        double previous = previousBodyLevelSets[node];
                        double current = lsm.LevelSetsBody[node];
                        string comparison = (previous == current) ? "SAME" : "DIFFERENT";
                        writer.WriteLine(
                            $"{node}: {comparison} - previous body level set = {previous} - current body level set = {current}");
                    }
                    writer.WriteLine();

                    writer.WriteLine("Old tip nodes, tip level set:");
                    foreach (var node in lsm.CrackTipNodesOld[lsm.CrackTipEnrichments])
                    {
                        double previous = previousTipLevelSets[node];
                        double current = lsm.LevelSetsTip[node];
                        string comparison = (previous == current) ? "SAME" : "DIFFERENT";
                        writer.WriteLine(
                            $"{node}: {comparison} - previous tip level set = {previous} - current tip level set = {current}");
                    }
                    writer.WriteLine();

                    writer.WriteLine("New tip nodes, tip level set:");
                    foreach (var node in lsm.CrackTipNodesNew[lsm.CrackTipEnrichments])
                    {
                        double previous = previousTipLevelSets[node];
                        double current = lsm.LevelSetsTip[node];
                        string comparison = (previous == current) ? "SAME" : "DIFFERENT";
                        writer.WriteLine(
                            $"{node}: {comparison} - previous tip level set = {previous} - current tip level set = {current}");
                    }
                    writer.WriteLine();
                }
            }
            previousBodyLevelSets = new Dictionary<XNode2D, double>(lsm.LevelSetsBody);
            previousTipLevelSets = new Dictionary<XNode2D, double>(lsm.LevelSetsTip);
            ++iteration;
        }
    }
}
