using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Logging.VTK
{
    /// <summary>
    /// Creates displacement, strain and stress observers that log the corresponding data to VTK output files. Then the .vtk
    /// files can be visualized in paraview.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class VtkLogFactory : ILogFactory
    {
        private readonly string directory;
        private readonly bool logDisplacements;
        private readonly bool logStrains;
        private readonly bool logStresses;
        private readonly Model model;

        public VtkLogFactory(Model model, string directory, bool logDisplacements, bool logStrains, bool logStresses)
        {
            this.model = model;
            this.directory = directory;

            if ((!logDisplacements) && (!logStrains) && (!logStresses)) throw new ArgumentException("You must log something");
            this.logDisplacements = logDisplacements;
            this.logStrains = logStrains;
            this.logStresses = logStresses;
        }

        public string DisplacementsFilename { get; set; } = "displacements";
        public string StrainsFilename { get; set; } = "strains";
        public string StressesFilename { get; set; } = "stresses";

        public IAnalyzerLog[] CreateLogs()
        {
            var logs = new List<IAnalyzerLog>(3);
            var mesh = new VtkMesh2D(model);
            if (logDisplacements) logs.Add(new Displacement2DLog(directory, DisplacementsFilename, model, mesh));
            //if (logStrains) logs.Add(new Strain2DLog(directory, DisplacementsFilename, model, mesh));
            //if (logStresses) logs.Add(new Stress2DLog(directory, DisplacementsFilename, model, mesh));
            return logs.ToArray();
        }
    }
}
