using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Output.VTK;

//TODO: Decide between 1) push observer, 2) pull observer with the observable injected in observer.Observe(observarble) and  
//      possibly generics for specific concrete observables or 3) pull observer with the observable injected during construction,
//      in which case, there is no need for generics but construction is harder and possibly needs builders.
namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.Logging
{
    class LevelSetLogger
    {
        private readonly TrackingExteriorCrackLSM lsm;
        private readonly Model2D model;
        private readonly string outputDirectory;
        private int iteration;

        public LevelSetLogger(Model2D model, TrackingExteriorCrackLSM lsm, string outputDirectory)
        {
            this.model = model;
            this.lsm = lsm;
            this.outputDirectory = outputDirectory;
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
            var crackWriter = new PolylineWriter();
            crackWriter.InitializeFile($"{outputDirectory}\\crack_{iteration}", true);
            crackWriter.WritePolyline(lsm.CrackPath);
            crackWriter.CloseCurrentFile();

            // Log the level sets
            var lsWriter = new VTKWriter(model);
            lsWriter.InitializeFile($"{outputDirectory}\\level_sets_{iteration}", true);

            IReadOnlyDictionary<XNode2D, double> levelSetsBody = lsm.LevelSetsBody;
            IReadOnlyDictionary<XNode2D, double> levelSetsTip = lsm.LevelSetsTip;
            var bodyArray = new double[levelSetsBody.Count];
            var tipArray = new double[levelSetsTip.Count];
            for (int i = 0; i < model.Nodes.Count; ++i)
            {
                XNode2D node = model.Nodes[i];
                bodyArray[i] = levelSetsBody[node];
                tipArray[i] = levelSetsTip[node];
            }

            lsWriter.WriteScalarField("level_set_crack_body", bodyArray);
            lsWriter.WriteScalarField("level_set_crack_tip", tipArray);

            lsWriter.CloseCurrentFile();
            ++iteration;
        }
    }
}
