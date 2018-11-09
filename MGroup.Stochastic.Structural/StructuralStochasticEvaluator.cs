using MGroup.Stochastic.Interfaces;
using MGroup.Stochastic.Structural.StochasticRealizers;
using System;
using System.Collections.Generic;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Analyzers;

namespace MGroup.Stochastic.Structural
{
    public class StructuralStochasticEvaluator : ISystemRealizer, ISystemResponseEvaluator
    {
        public double YoungModulus { get; }
        public IStochasticDomainMapper DomainMapper;
        //public RandomVariable StochasticRealization { get; }
        public KarhunenLoeveCoefficientsProvider StochasticRealization { get; }
        //public ModelBuilder ModelBuilder { get; }
        public GiannisModelBuilder ModelBuilder { get; }
        private Model currentModel;
        int karLoeveTerms = 2;
        double[] domainBounds = new double[2] { 0, 1.0 };
        double sigmaSquare = 0.1;
        double meanValue = 1;
        int partition = 11;
        double correlationLength = 10.0;
        bool isGaussian = true;
        int PCorder = 1;
        bool midpointMethod = true;
        int mcsamples = 5;

        //public StructuralStochasticEvaluator(double youngModulus, IStochasticDomainMapper domainMapper)
        //{
        //    YoungModulus = youngModulus;
        //    DomainMapper = domainMapper;
        //    ModelBuilder = new ModelBuilder(); 
        //    StochasticRealization = new RandomVariable(youngModulus, domainMapper);
        //}

        public StructuralStochasticEvaluator(double youngModulus, IStochasticDomainMapper domainMapper)
        {
            YoungModulus = youngModulus;
            DomainMapper = domainMapper;
            ModelBuilder = new GiannisModelBuilder();
            StochasticRealization = new KarhunenLoeveCoefficientsProvider(partition, youngModulus, midpointMethod,
                isGaussian, karLoeveTerms, domainBounds, sigmaSquare, correlationLength);
        }

        public void Realize(int iteration)
        {
            currentModel = ModelBuilder.GetModel(StochasticRealization, DomainMapper, iteration);
        }



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
            return new[] { linearSystems[0].RHS[0] };
        }

    }
}
