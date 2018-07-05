using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.VonMisesStress;
using System;

//TODO: this should take as input the PreprocessorModel and get the problem type and materials from there.
namespace ISAAR.MSolve.Preprocessor.UI
{
    /// <summary>
    /// Utility class to facilitate the definition the field (displacements, stresses, etc) output requests. The requested  
    /// fields will be gathered during the simulation and exported to appropriate formats. The user may process the results 
    /// once the simulation ends.
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

        internal ILogFactory CreateLogFactory(PreprocessorModel model)
        {
            var logFactory = new VtkLogFactory(model.CoreModel, outputDirectory)
            {
                LogDisplacements = Displacements,
                LogStrains = Strains,
                LogStresses = Stresses,
            };

            if (StressesVonMises)
            {
                if (model.Dimensions == PreprocessorModel.ProblemDimensions.TwoDimensionalPlaneStress)
                {
                    logFactory.VonMisesStressCalculator = new PlaneStressVonMises();
                }
                else if (model.Dimensions == PreprocessorModel.ProblemDimensions.TwoDimensionalPlaneStrain)
                {
                    logFactory.VonMisesStressCalculator = new ElasticPlaneStrainVonMises(model.PlainStrainMaterial);
                }
                else if (model.Dimensions == PreprocessorModel.ProblemDimensions.ThreeDimensional)
                {
                    throw new NotImplementedException();
                }
            }

            return logFactory;
        }
    }
}