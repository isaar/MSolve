using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Postprocessing;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Logging
{
    /// <summary>
    /// Observer that logs the displacement field at each analysis step to .vtk output files.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Displacement2DLog : IAnalyzerLog
    {
        private readonly DisplacementField2D displacementField;
        private readonly VtkMesh2D mesh;
        private readonly Model model;
        private readonly string pathNoExtension;

        //TODO: This should be controlled by the Analyzer and passed in through IAnalyzerLog.StoreResults().
        private int iteration;

        /// <summary>
        /// Instantiates a new log that writes the displacement field of each iteration to output files, which can then be 
        /// processed in Paraview.
        /// </summary>
        /// <param name="directory">The full path of the folder whre the output files will be written, e.g. 
        ///     C:\\Users\\MyUser\\Desktop\\Paraview.</param>
        /// <param name="filename">The name of the displacement file(s) without extensions, e.g. displacements. Keep in mind
        ///     that the actual files will be suffixed with the iteration number, e.g. 
        ///     C:\\Users\\MyUser\\Desktop\\Paraview\\displacements_0.vtk, 
        ///     C:\\Users\\MyUser\\Desktop\\Paraview\\displacements_1.vtk, etc.</param>
        /// <param name="model">A collection of nodes, elements and other entities.</param>
        /// <param name="mesh">The same mesh that is contained in <paramref name="model"/>, but expressed in VTK objects. This 
        ///     object should be shared across other VTK logs to reduce memory consumption.</param>
        public Displacement2DLog(string directory, string filename, Model model, VtkMesh2D mesh)
        {
            this.pathNoExtension = directory + "\\" + filename;
            this.model = model;
            this.mesh = mesh;
            this.displacementField = new DisplacementField2D(model);
            iteration = 0;
        }

        public void StoreResults(DateTime startTime, DateTime endTime, IVector solution)
        {
            // Find displacements of each VTK point
            var displacements = new double[mesh.Points.Count][];
            displacementField.FindNodalDisplacements(solution);
            foreach (Node node in model.Nodes)
            {
                int pointID = mesh.Nodes2Points[node].ID;
                displacements[pointID] = displacementField[node];
            }

            // Write them to file
            string path = pathNoExtension + $"_{iteration}.vtk";
            using (var writer = new VtkFileWriter(path))
            {
                writer.WriteMesh(mesh.Points, mesh.Cells);
                writer.WriteVector2DField("displacements", displacements);
            }

            ++iteration;
        }
    }
}
