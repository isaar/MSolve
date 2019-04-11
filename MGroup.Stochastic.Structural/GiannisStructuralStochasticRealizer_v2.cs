using ISAAR.MSolve.FEM.Entities;
using MGroup.Stochastic.Interfaces;
using MGroup.Stochastic.Structural.StochasticRealizers;

namespace MGroup.Stochastic.Structural
{
    public class GiannisStructuralStochasticRealizer_v2 : ISystemRealizer
    {
        public double YoungModulus { get; }
        private readonly IStochasticDomainMapper DomainMapper;
        public KarhunenLoeveCoefficientsProvider StochasticRealization { get; }
        public GiannisModelBuilder_v2 ModelBuilder { get; }
        private Model_v2 currentModel;
        int karLoeveTerms = 2;
        double[] domainBounds = new double[2] { 0, 1.0 };
        double sigmaSquare = 0.01;
        double meanValue = 1;
        int partition = 11;
        double correlationLength = 1.0;
        bool isGaussian = true;
        int PCorder = 1;
        bool midpointMethod = true;
        int mcsamples = 5;


        /// <summary>Initializes a new instance of the <see cref="GiannisStructuralStochasticRealizer"/> class.</summary>
        /// <param name="youngModulus">The young modulus.</param>
        /// <param name="domainMapper">The domain mapper.</param>
        public GiannisStructuralStochasticRealizer_v2(double youngModulus, IStochasticDomainMapper domainMapper)
        {
            YoungModulus = youngModulus;
            meanValue = youngModulus;
            DomainMapper= domainMapper;
            ModelBuilder = new GiannisModelBuilder_v2();
            StochasticRealization = new KarhunenLoeveCoefficientsProvider(partition, meanValue, midpointMethod, 
            isGaussian, karLoeveTerms, domainBounds, sigmaSquare, correlationLength);
        }

        /// <summary>Realizes the specified iteration.</summary>
        /// <param name="iteration">The iteration.</param>
        public void Realize(int iteration)
        {
            currentModel = ModelBuilder.GetModel(StochasticRealization, DomainMapper, iteration);
        }
    }
}
