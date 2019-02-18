using System;
using System.Collections.Generic;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Stochastic;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;
using Troschuetz.Random.Distributions.Continuous;

namespace ISAAR.MSolve.Analyzers
{
    public class MonteCarloAnalyzerStiffnessMatrixFactoryWithStochasticMaterial_v2 : IParentAnalyzer
    {
        private int currentSimulation = -1;
        private readonly int expansionOrder;
        private readonly int simulations;
        private readonly IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems;
        //private readonly IDictionary<int, IMatrix<double>> matrices;
        private readonly IDictionary<int, IMatrix>[] matrices;
        private readonly Model_v2 model;
        private readonly IAnalyzerProvider_v2 provider;
        private readonly double[][] randomNumbers;
        private readonly IStochasticMaterialCoefficientsProvider coefficientsProvider;

        public MonteCarloAnalyzerStiffnessMatrixFactoryWithStochasticMaterial_v2(Model_v2 model, IAnalyzerProvider_v2 provider, 
            IChildAnalyzer embeddedAnalyzer, ISolver_v2 solver, 
            GaussianFileStochasticCoefficientsProvider coefficientsProvider, int expansionOrder, int simulations)
        {
            this.ChildAnalyzer = embeddedAnalyzer;
            this.provider = provider;
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.expansionOrder = expansionOrder;
            this.simulations = simulations;
            this.ChildAnalyzer.ParentAnalyzer = this;
            //this.matrices = new Dictionary<int, IMatrix<double>>(subdomains.Count);
            this.matrices = new Dictionary<int, IMatrix>[expansionOrder + 1];
            this.randomNumbers = new double[simulations][];
            this.coefficientsProvider = coefficientsProvider;
            //this.stochasticDomain = stochasticDomain;

            NormalDistribution n = new NormalDistribution();
            n.Mu = 0;
            n.Sigma = 1;
            string[] randoms = new string[simulations];
            for (int i = 0; i < simulations; i++)
            {
                randomNumbers[i] = new double[expansionOrder];
                for (int j = 0; j < expansionOrder; j++)
                    randomNumbers[i][j] = n.NextDouble();
                randoms[i] = randomNumbers[i][0].ToString();
            }
            //File.WriteAllLines(String.Format(@"randoms.txt", expansionOrder), randoms);
        }

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog_v2[]> Logs { get; } = new Dictionary<int, IAnalyzerLog_v2[]>();

        public IChildAnalyzer ChildAnalyzer { get; set; }

        public void BuildMatrices()
        {
            if (ChildAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");
            if (currentSimulation < 0) return;

            provider.Reset();
            ChildAnalyzer.BuildMatrices();
        }

        public void Initialize(bool isFirstAnalysis)
        {
            if (ChildAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");

            for (int i = 0; i < simulations; i++)
            {
                currentSimulation = i;
                coefficientsProvider.RandomVariables = randomNumbers[currentSimulation];
                BuildMatrices();
            }
        }

        public void Solve()
        {
            throw new InvalidOperationException("Monte Carlo stifness factory cannot solve.");
        }

        #endregion
    }
}
