using MGroup.Stochastic.Interfaces;
using MGroup.Stochastic.Structural.StochasticRealizers;
using System;
using System.Collections.Generic;
using Accord;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace MGroup.Stochastic.Structural
{
    public class StructuralStochasticEvaluator : ISystemRealizer, ISystemResponseEvaluator
    {
        public double YoungModulus { get; }
        public IStochasticDomainMapper DomainMapper;
        //public RandomVariable StochasticRealization { get; }
        public KarhunenLoeveCoefficientsProvider StochasticRealization { get; }
        //public SpectralRepresentation1DRandomFieldGenerator StochasticRealization { get; }
        //public ModelBuilder ModelBuilder { get; }
        public GiannisModelBuilder ModelBuilder { get; }
        private Model currentModel;
        int karLoeveTerms = 4;
        double[] domainBounds = new double[2] { 0, 1 };
        double sigmaSquare = .01;
        double meanValue = 1;
        int partition = 21;
        double correlationLength = 1.0;
        bool isGaussian = true;
        int PCorder = 1;
        bool midpointMethod = true;

        //public StructuralStochasticEvaluator(double youngModulus, IStochasticDomainMapper domainMapper)
        //{
        //    YoungModulus = youngModulus;
        //    DomainMapper = domainMapper;
        //    ModelBuilder = new ModelBuilder(); 
        //    StochasticRealization = new RandomVariable(youngModulus, domainMapper);
        //}

        /// <summary>Initializes a new instance of the <see cref="StructuralStochasticEvaluator"/> class.</summary>
        /// <param name="youngModulus">The young modulus.</param>
        /// <param name="domainMapper">The domain mapper.</param>
        public StructuralStochasticEvaluator(double youngModulus, IStochasticDomainMapper domainMapper)
        {
            YoungModulus = youngModulus;
            DomainMapper = domainMapper;
            ModelBuilder = new GiannisModelBuilder();
            StochasticRealization = new KarhunenLoeveCoefficientsProvider(partition, youngModulus, midpointMethod,
                isGaussian, karLoeveTerms, domainBounds, sigmaSquare, correlationLength);
            //StochasticRealization = new SpectralRepresentation1DRandomFieldGenerator(10, 0.1, youngModulus, .05, 0.1, 256);
        }

        /// <summary>Realizes the specified iteration.</summary>
        /// <param name="iteration">The iteration.</param>
        public void Realize(int iteration)
        {
            currentModel = ModelBuilder.GetModel(StochasticRealization, DomainMapper, iteration);
        }



        /// <summary>Evaluates the specified iteration.</summary>
        /// <param name="iteration">The iteration.</param>
        /// <returns></returns>
        public double[] Evaluate(int iteration)
        {
            var linearSystems = new Dictionary<int, ILinearSystem>
            {
                { 0, new SkylineLinearSystem(0, currentModel.SubdomainsDictionary[0].Forces) }
            };
            var solver = new SolverSkyline(linearSystems[0]);
            var provider = new ProblemStructural(currentModel, linearSystems);
            var childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            var parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            //return new[] { linearSystems[0].RHS[0] };
            return new[] { linearSystems[0].Solution[56], linearSystems[0].Solution[58] };
        }

    }
}
