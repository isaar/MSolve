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
        public RandomVariable StochasticRealization { get; }
        public ModelBuilder ModelBuilder { get; }
        private Model currentModel;

        public StructuralStochasticEvaluator(double youngModulus, IStochasticDomainMapper domainMapper)
        {
            YoungModulus = youngModulus;
            DomainMapper = domainMapper;
            ModelBuilder = new ModelBuilder(); 
            StochasticRealization = new RandomVariable(youngModulus, domainMapper);
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
