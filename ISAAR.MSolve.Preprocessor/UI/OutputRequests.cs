using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.VonMisesStress;

namespace ISAAR.MSolve.Preprocessor.UI
{
    public class OutputRequests
    {
        private readonly string outputDirectory;

        public OutputRequests(string outputDirectory)
        {
            this.outputDirectory = outputDirectory;
        }

        public bool Displacements { get; set; } = true;
        public bool Strains { get; set; } = true;
        public bool Stresses { get; set; } = true;
        public ElasticMaterial2D Material { get; set; } //TODO: this should be taken by the model's representation

        internal ILogFactory CreateLogFactory(Model model)
        {
            var logFactory = new VtkLogFactory(model, outputDirectory)
            {
                LogDisplacements = Displacements,
                LogStrains = Strains,
                LogStresses = Stresses,
            };

            if (Stresses)
            {
                if (Material.StressState == StressState2D.PlaneStress)
                {
                    logFactory.VonMisesStressCalculator = new PlaneStressVonMises();
                }
                else
                {
                    logFactory.VonMisesStressCalculator = new ElasticPlaneStrainVonMises(Material);
                }
            }

            return logFactory;
        }
    }
}