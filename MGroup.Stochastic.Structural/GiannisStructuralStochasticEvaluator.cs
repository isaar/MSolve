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
    public class GiannisStructuralStochasticEvaluator : ISystemRealizer
    {
        public double YoungModulus { get; }
        public KarhunenLoeveCoefficientsProvider StochasticRealization { get; }
        public ModelBuilder ModelBuilder { get; }
        private Model currentModel;
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

        
        public GiannisStructuralStochasticEvaluator(double youngModulus)
        {
            YoungModulus = youngModulus;
            ModelBuilder = new ModelBuilder();
            StochasticRealization = new KarhunenLoeveCoefficientsProvider( mcsamples, partition, meanValue, midpointMethod, 
            isGaussian, karLoeveTerms, domainBounds, sigmaSquare, correlationLength);

            double[] xCoordinates = StochasticRealization.KarhunenLoeveFredholmWithFEM(KarLoeveTerms, domainBounds, sigmaSquare, partition, correlationLength).Item1;
            double[] lambda = StochasticRealization.KarhunenLoeveFredholmWithFEM(KarLoeveTerms, domainBounds, sigmaSquare, partition, correlationLength).Item2;
            double[,] Eigenvectors = StochasticRealization.KarhunenLoeveFredholmWithFEM(KarLoeveTerms, domainBounds, sigmaSquare, partition, correlationLength).Item3;
            double[,] fieldRealizations = StochasticRealization.KarhunenLoeveFredholm1DSampleGenerator(MCsamples, lambda, Eigenvectors, meanValue, midpointMethod, isGaussian);
        }

        public void Realize(int iteration)
        {
            currentModel = ModelBuilder.GetModel(StochasticRealization, iteration);
        }

        public double[] Evaluate(int iteration)
        {
            var linearSystems = new Dictionary<int, ILinearSystem>
            {
                { 0, new SkylineLinearSystem(0, currentModel.SubdomainsDictionary[0].Forces) }
            };
            var s = new SolverSkyline(linearSystems[0]);
            var p = new ProblemStructural(currentModel, linearSystems);
            var a = new LinearAnalyzer(s, linearSystems);
            var sa = new StaticAnalyzer(p, a, linearSystems);

            sa.Initialize();
            sa.Solve();
            return new[] { linearSystems[0].RHS[0] };
        }

    }
}
