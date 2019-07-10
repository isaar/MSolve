using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.XFEM.Output.VTK;
using ISAAR.MSolve.Logging.VTK;

//TODO: Decide between 1) push observer, 2) pull observer with the observable injected in observer.Observe(observarble) and  
//      possibly generics for specific concrete observables or 3) pull observer with the observable injected during construction,
//      in which case, there is no need for generics but construction is harder and possibly needs builders.
namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.Logging
{
    public class LevelSetLogger
    {
        private readonly TrackingExteriorCrackLSM lsm;
        private readonly XModel model;
        private readonly string outputDirectory;
        private readonly VtkMesh<XNode> vtkMesh;
        private int iteration;

        public LevelSetLogger(XModel model, TrackingExteriorCrackLSM lsm, string outputDirectory)
        {
            this.model = model;
            this.lsm = lsm;
            this.outputDirectory = outputDirectory;
            this.vtkMesh = new VtkMesh<XNode>(model.Nodes, model.Elements);
        }

        // This could be handled better by having LSM storing the crack path (which is all around useful). However having an
        // InitialLog() method also provides some versatility.
        public void InitialLog()
        {
            iteration = 0;
            Log();
        }

        public void Log() 
        {
            // Log the crack path
            using (var crackWriter = new VtkPolylineWriter($"{outputDirectory}\\crack_{iteration}.vtk"))
            {
                crackWriter.WritePolyline(lsm.CrackPath);
            }

            // Log the level sets
            using (var lsWriter = new VtkFileWriter($"{outputDirectory}\\level_sets_{iteration}.vtk"))
            {
                // The mesh should be provided by a dedicated IVtkMesh. Its implementations will be continuous/discontinuous, 
                // FEM/XFEM. These classes should then replace the relevant code in DiscontinuousMeshVTKWriter.cs and Logging project.
                // Also the Dictionary<Element, VtkCell> should be avoided. Just store the VtkCells in the same order as the IList<Element> of the model.
                // Ideally VtkPoint can be completely replaced by INode and perhaps VtkCell can be replaced by IElement or ICell<TNode>.
                lsWriter.WriteMesh(vtkMesh);
                IReadOnlyDictionary<XNode, double> levelSetsBody = lsm.LevelSetsBody;
                IReadOnlyDictionary<XNode, double> levelSetsTip = lsm.LevelSetsTip;
                var bodyArray = new double[levelSetsBody.Count];
                var tipArray = new double[levelSetsTip.Count];
                for (int i = 0; i < model.Nodes.Count; ++i)
                {
                    XNode node = model.Nodes[i];
                    bodyArray[i] = levelSetsBody[node];
                    tipArray[i] = levelSetsTip[node];
                }

                lsWriter.WriteScalarField("level_set_crack_body", bodyArray);
                lsWriter.WriteScalarField("level_set_crack_tip", tipArray);
            }
            ++iteration;
        }
    }
}
