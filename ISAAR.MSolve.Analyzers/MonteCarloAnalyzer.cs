using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using Troschuetz.Random.Distributions.Continuous;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Stochastic;

namespace ISAAR.MSolve.Analyzers
{
    public class MonteCarloAnalyzer : IAnalyzer
    {
        private int currentSimulation = -1;
        private readonly int expansionOrder;
        private readonly int simulations;
        private readonly IDictionary<int, ILinearSystem> subdomains;
        //private readonly IDictionary<int, IMatrix2D<double>> matrices;
        private readonly IDictionary<int, IMatrix2D>[] matrices;
        private readonly Model model;
        private readonly Dictionary<int, IAnalyzerLog[]> logs = new Dictionary<int, IAnalyzerLog[]>();
        private readonly IAnalyzerProvider provider;
        private readonly double[][] randomNumbers;
        //private readonly double[] stochasticDomain;
        private IAnalyzer childAnalyzer;
        private IAnalyzer parentAnalyzer = null;
        private readonly GaussianFileStochasticCoefficientsProvider coefficientsProvider;

        //public MonteCarloAnalyzer(Model model, IAnalyzerProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ISolverSubdomain> subdomains, double[] stochasticDomain, int expansionOrder, int simulations)
        public MonteCarloAnalyzer(Model model, IAnalyzerProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ILinearSystem> subdomains, GaussianFileStochasticCoefficientsProvider coefficientsProvider, int expansionOrder, int simulations)
        {
            this.childAnalyzer = embeddedAnalyzer;
            this.provider = provider;
            this.model = model;
            this.subdomains = subdomains;
            this.expansionOrder = expansionOrder;
            this.simulations = simulations;
            this.childAnalyzer.ParentAnalyzer = this;
            //this.matrices = new Dictionary<int, IMatrix2D<double>>(subdomains.Count);
            this.matrices = new Dictionary<int, IMatrix2D>[expansionOrder + 1];
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
            File.WriteAllLines(String.Format(@"randoms.txt", expansionOrder), randoms);
        }

        private void InitializeCoefficientsProvider()
        {
            //coefficientsProvider = new FileStochasticCoefficientsProvider("Lognormal.csv", stochasticDomain);
            //coefficientsProvider = new FileStochasticCoefficientsProvider("Gaussian.txt", 50, '\t', stochasticDomain);
            foreach (var subdomain in model.Subdomains)
                foreach (var e in subdomain.ElementsDictionary.Values.Where(e => e.ElementType is IStochasticFiniteElement))
                    ((IStochasticFiniteElement)e.ElementType).CoefficientsProvider = coefficientsProvider;
        }

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog[]> Logs { get { return logs; } }
        public IAnalyzer ParentAnalyzer
        {
            get { return parentAnalyzer; }
            set { parentAnalyzer = value; }
        }

        public IAnalyzer ChildAnalyzer
        {
            get { return childAnalyzer; }
            set { childAnalyzer = value; }
        }

        public void BuildMatrices()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");
            if (currentSimulation < 0)
                InitializeCoefficientsProvider();

            provider.Reset();
            coefficientsProvider.CurrentOrder = -1;
            childAnalyzer.BuildMatrices();

            matrices[0] = new Dictionary<int, IMatrix2D>(subdomains.Count);
            foreach (var subdomain in subdomains.Values)
            {
                SkylineMatrix2D k = (SkylineMatrix2D)subdomain.Matrix;
                matrices[0].Add(subdomain.ID, (SkylineMatrix2D)k.Clone());
            }
            for (int i = 0; i < expansionOrder; i++)
            {
                provider.Reset();
                coefficientsProvider.CurrentOrder = i;
                childAnalyzer.BuildMatrices();

                matrices[i + 1] = new Dictionary<int, IMatrix2D>(subdomains.Count);
                foreach (var subdomain in subdomains.Values)
                {
                    SkylineMatrix2D k = (SkylineMatrix2D)subdomain.Matrix;
                    matrices[i + 1].Add(subdomain.ID, (SkylineMatrix2D)k.Clone());
                }
            }

            //matrices.Clear();
            //provider.Reset();
            //coefficientsProvider.CurrentOrder = -1;
            //childAnalyzer.BuildMatrices();

            //foreach (var subdomain in subdomains.Values)
            //{
            //    SkylineMatrix2D<double> k = (SkylineMatrix2D<double>)subdomain.Matrix;
            //    matrices.Add(subdomain.ID, (SkylineMatrix2D<double>)k.Clone());
            //}
            //for (int i = 0; i < expansionOrder; i++)
            //{
            //    provider.Reset();
            //    coefficientsProvider.CurrentOrder = i;
            //    childAnalyzer.BuildMatrices();

            //    foreach (var subdomain in subdomains.Values)
            //    {
            //        SkylineMatrix2D<double> k = (SkylineMatrix2D<double>)subdomain.Matrix;
            //        matrices[subdomain.ID].LinearCombination(new double[] { randomNumbers[currentSimulation][i] }, new IMatrix2D<double>[] { k });
            //    }
            //}
            //foreach (var subdomain in subdomains.Values)
            //    subdomain.Matrix = matrices[subdomain.ID];
        }

        private void ComposeStochasticMatrixFromMatrices()
        {
            var currentRandomNumbers = randomNumbers[currentSimulation];
            var coefficients = new double[] { 1 }.Concat(currentRandomNumbers).ToList<double>();
            var matricesPerSubdomain = new Dictionary<int, IMatrix2D[]>();
            foreach (var subdomain in subdomains.Values)
            {
                int id = subdomain.ID;
                var tempMatrices = new IMatrix2D[expansionOrder + 1];
                for (int i = 0; i <= expansionOrder; i++)
                    tempMatrices[i] = matrices[i][id];
                matricesPerSubdomain.Add(id, tempMatrices);
            }

            foreach (var subdomain in subdomains.Values)
                subdomain.Matrix = (SkylineMatrix2D)((SkylineMatrix2D)matrices[0][subdomain.ID]).Clone();
            foreach (var subdomain in subdomains.Values)
                ((ILinearlyCombinable)subdomain.Matrix).LinearCombination(coefficients, matricesPerSubdomain[subdomain.ID]);

            
            //for (int i = 0; i < expansionOrder; i++)
            //{
            //    foreach (var subdomain in subdomains.Values)
            //    {
            //        SkylineMatrix2D<double> k = (SkylineMatrix2D<double>)matrices[i + 1][subdomain.ID];
            //        subdomain.Matrix.LinearCombination(new double[] { randomNumbers[currentSimulation][i] }, new IMatrix2D<double>[] { k });
            //    }
            //}
        }

        public void Initialize()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");
        }

        public void Solve()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");

            string[] values = new string[simulations];
            for (int i = 0; i < simulations; i++)
            {
                currentSimulation = i;
                ComposeStochasticMatrixFromMatrices();
                childAnalyzer.Initialize();
                childAnalyzer.Solve();
                values[i] = subdomains[1].Solution[28].ToString();
            }

            File.WriteAllLines(String.Format(@"montecarlo{0}.txt", expansionOrder), values);
        }

        #endregion
    }
}
