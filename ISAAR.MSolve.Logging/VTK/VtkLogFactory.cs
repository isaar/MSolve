using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Logging.VTK
{
    /// <summary>
    /// Creates <see cref="VtkLog2D"/> observers that log the corresponding data to VTK output files. Then the .vtk
    /// files can be visualized in paraview.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class VtkLogFactory : ILogFactory
    {
        private readonly string directory;
        private readonly Model model;

        public VtkLogFactory(Model model, string directory)
        {
            this.model = model;
            this.directory = directory;
        }

        public bool LogDisplacements { get; set; } = true;
        public bool LogStrains { get; set; } = true;
        public bool LogStresses { get; set; } = true;
        public string Filename { get; set; } = "field_output";

        public IAnalyzerLog[] CreateLogs()
        {
            var logs = new List<IAnalyzerLog>(3);
            var mesh = new VtkMesh2D(model);
            logs.Add(new VtkLog2D(directory, Filename, model, mesh, LogDisplacements, LogStrains, LogStresses));
            return logs.ToArray();
        }
    }
}
