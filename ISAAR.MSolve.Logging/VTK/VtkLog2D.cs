using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Postprocessing;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Materials.VonMisesStress;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Logging.VTK
{
    /// <summary>
    /// Observer that logs the displacement, strain and stress field at each analysis step to .vtk output files (1 file per 
    /// analysis step).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class VtkLog2D : IAnalyzerLog
    {
        private readonly VtkMesh2D mesh;
        private readonly Model model;
        private readonly string pathNoExtension;
        private readonly bool logDisplacements, logStrains, logStresses;
        private readonly IVonMisesStress2D vonMisesStressCalculator;

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
        /// <param name="logDisplacements">If true, the displacement field will also be written to the output file.</param>
        /// <param name="logStrains">If true, the strain field will also be written to the output file.</param>
        /// <param name="logStresses">If true, the stress field will also be written to the output file.</param>
        /// <param name="vonMisesStressCalculator">The strategy used for calculating von Mises equivalent stress. If null is 
        ///     passed, then von Mises equivalent stresses will not be written to the output file.</param>
        public VtkLog2D(string directory, string filename, Model model, VtkMesh2D mesh, 
            bool logDisplacements, bool logStrains, bool logStresses, IVonMisesStress2D vonMisesStressCalculator)
        {
            this.pathNoExtension = directory + "\\" + filename;
            this.model = model;
            this.mesh = mesh;
            this.logDisplacements = logDisplacements;
            this.logStrains = logStrains;
            this.logStresses = logStresses;
            this.vonMisesStressCalculator = vonMisesStressCalculator;
            iteration = 0;
        }

        public void StoreResults(DateTime startTime, DateTime endTime, IVector solution)
        {
            string path = pathNoExtension + $"_{iteration}.vtk";
            using (var writer = new VtkFileWriter(path))
            {
                writer.WriteMesh(mesh.Points, mesh.Cells);

                //TODO: this should be abstracted, especially the nodes.
                IList<Node> nodes = model.Nodes;
                int numPoints = mesh.Points.Count;
                if (nodes.Count != numPoints)
                {
                    throw new NotImplementedException("For now the nodes and VTK mesh vertices must be 1-1");
                }
                var pointIDs = new int[nodes.Count];
                for (int i = 0; i < nodes.Count; ++i) pointIDs[i] = mesh.Nodes2Points[nodes[i]].ID;

                // Find the displacements of each VTK point and write them to the output file
                if (logDisplacements)
                {
                    var displacementField = new DisplacementField2D(model);
                    displacementField.FindNodalDisplacements(solution);
                    var displacements = new double[numPoints][]; //TODO: this conversion between data structures should be avoided. Redundant memory and computations.
                    for (int i = 0; i < numPoints; ++i) displacements[pointIDs[i]] = displacementField[nodes[i]];
                    writer.WriteVector2DField("displacements", displacements);
                }
                
                if (logStrains || logStresses || (vonMisesStressCalculator != null))
                {
                    var tensorsField = new StrainStressField2D(model);
                    tensorsField.CalculateNodalTensors(solution);

                    if (logStrains)
                    {
                        var strains = new double[numPoints][]; //TODO: this conversion between data structures should be avoided. Redundant memory and computations.
                        for (int i = 0; i < numPoints; ++i) strains[pointIDs[i]] = tensorsField.GetStrainsOfNode(nodes[i]);
                        writer.WriteTensor2DField("strain", strains);
                    }

                    if (logStresses)
                    {
                        var stresses = new double[numPoints][];
                        for (int i = 0; i < numPoints; ++i) stresses[pointIDs[i]] = tensorsField.GetStressesOfNode(nodes[i]);
                        writer.WriteTensor2DField("stress", stresses);

                    }

                    //TODO: von Mises stress should probably be computed per Gauss point, then follow the usual extrapolation 
                    //      and averaging process
                    if (vonMisesStressCalculator != null)
                    {
                        var vonMises = new double[numPoints];
                        for (int i = 0; i < numPoints; ++i)
                        {
                            vonMises[pointIDs[i]] = vonMisesStressCalculator.Calculate(tensorsField.GetStrainsOfNode(nodes[i]), 
                                tensorsField.GetStressesOfNode(nodes[i]));
                        }
                        writer.WriteScalarField("stress_vonMises", vonMises);
                    }
                }
            }

            ++iteration;
        }
    }
}
