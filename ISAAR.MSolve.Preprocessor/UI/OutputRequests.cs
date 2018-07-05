using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.VonMisesStress;

//TODO: this should take as input the PreprocessorModel and get the problem type and materials from there.
namespace ISAAR.MSolve.Preprocessor.UI
{
    /// <summary>
    /// Defines the field (displacements, stresses, etc) output requests. The requested fields will be gathered during the 
    /// simulation and exported to appropriate formats. The user may process the results once the simulation ends.
    /// </summary>
    public class OutputRequests
    {
        private readonly string outputDirectory;

        /// <summary>
        /// Instantiates a new <see cref="OutputRequests"/>.
        /// </summary>
        /// <param name="outputDirectory">The absolute path of the directory where output files will be written to.</param>
        public OutputRequests(string outputDirectory)
        {
            this.outputDirectory = outputDirectory;
        }

        public bool Displacements { get; set; } = true;

        /// <summary>
        /// Engineering strain tensor along the axes of the model. The values at nodes where more than one elements intersect 
        /// will be averaged. 
        /// </summary>
        public bool Strains { get; set; } = false;

        /// <summary>
        /// Cauchy stress tensor along the axes of the model. The values at nodes where more than one elements intersect will be 
        /// averaged. 
        /// </summary>
        public bool Stresses { get; set; } = false;

        /// <summary>
        /// The equivalent von Mises stress that can be used to check the von Mises yield criterion.
        /// </summary>
        public bool StressesVonMises { get; set; } = false;

        /// <summary>
        /// Defines whether the problem is <see cref="StressState2D.PlaneStress"/> or <see cref="StressState2D.PlaneStrain"/>.
        /// Not applicable to 3D problems.
        /// </summary>
        public StressState2D StressState { get; set; } = StressState2D.PlaneStress; //TODO: this should be taken by the model's representation

        /// <summary>
        /// The material properties. Only applicable to <see cref="StressState2D.PlaneStrain"/> problems.
        /// </summary>
        public ElasticMaterial2D PlaneStrainMaterial { get; set; } = null; //TODO: this should be taken by the model's representation

        internal ILogFactory CreateLogFactory(Model model)
        {
            var logFactory = new VtkLogFactory(model, outputDirectory)
            {
                LogDisplacements = Displacements,
                LogStrains = Strains,
                LogStresses = Stresses,
            };

            if (StressesVonMises)
            {
                if (StressState == StressState2D.PlaneStress)
                {
                    logFactory.VonMisesStressCalculator = new PlaneStressVonMises();
                }
                else
                {
                    logFactory.VonMisesStressCalculator = new ElasticPlaneStrainVonMises(PlaneStrainMaterial);
                }
            }

            return logFactory;
        }
    }
}