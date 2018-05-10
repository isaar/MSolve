using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Output.VTK;

//TODO: Decide between 1) push observer, 2) pull observer with the observable injected in observer.Observe(observarble) and  
//      possibly generics for specific concrete observables or 3) pull observer with the observable injected during construction,
//      in which case, there is no need for generics but construction is harder and possibly needs builders.
namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    class LsmLogger
    {
        private readonly Model2D model;
        private readonly TrackingExteriorCrackLSM lsm;
        private readonly string outputDirectory;
        private int iteration;

        public LsmLogger(Model2D model, TrackingExteriorCrackLSM lsm, string outputDirectory)
        {
            this.model = model;
            this.lsm = lsm;
            this.outputDirectory = outputDirectory;
            this.iteration = 0;
        }

        public void Log() 
        {
            // TODO: Log the crack path


            // Log the level sets
            var writer = new VTKWriter(model);
            writer.InitializeFile($"{outputDirectory}\\level_sets_{iteration}", true);

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

            writer.WriteScalarField("level_set_crack_body", bodyArray);
            writer.WriteScalarField("level_set_crack_tip", tipArray);

            writer.CloseCurrentFile();
            ++iteration;
        }
    }
}
