using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Stochastic;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;
using Troschuetz.Random.Distributions.Continuous;

namespace ISAAR.MSolve.Analyzers
{
    public class MonteCarloAnalyzer_v2 : IParentAnalyzer
    {
        private int currentSimulation = -1;
        private readonly int expansionOrder;
        private readonly int simulations;
        private readonly IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems;
        private readonly IDictionary<int, IMatrix>[] matrices;
        private readonly Model_v2 model;
        private readonly Dictionary<int, IAnalyzerLog_v2[]> logs = new Dictionary<int, IAnalyzerLog_v2[]>();
        private readonly IAnalyzerProvider_v2 provider;
        private readonly double[][] randomNumbers;
        private readonly GaussianFileStochasticCoefficientsProvider coefficientsProvider;

        //public MonteCarloAnalyzer(Model model, IAnalyzerProvider_v2 provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ISolverSubdomain> subdomains, double[] stochasticDomain, int expansionOrder, int simulations)
        public MonteCarloAnalyzer_v2(Model_v2 model, IAnalyzerProvider_v2 provider, IChildAnalyzer embeddedAnalyzer, 
            ISolver_v2 solver, GaussianFileStochasticCoefficientsProvider coefficientsProvider, 
            int expansionOrder, int simulations)
        {
            this.ChildAnalyzer = embeddedAnalyzer;
            this.provider = provider;
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.expansionOrder = expansionOrder;
            this.simulations = simulations;
            this.ChildAnalyzer.ParentAnalyzer = this;
            //this.matrices = new Dictionary<int, IMatrix2D<double>>(subdomains.Count);
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
            File.WriteAllLines(String.Format(@"randoms.txt", expansionOrder), randoms);
        }

        private void InitializeCoefficientsProvider()
        {
            //coefficientsProvider = new FileStochasticCoefficientsProvider("Lognormal.csv", stochasticDomain);
            //coefficientsProvider = new FileStochasticCoefficientsProvider("Gaussian.txt", 50, '\t', stochasticDomain);
            foreach (var subdomain in model.Subdomains)
                foreach (var e in subdomain.Elements.Where(e => e.ElementType is IStochasticFiniteElement_v2))
                    ((IStochasticFiniteElement_v2)e.ElementType).CoefficientsProvider = coefficientsProvider;
        }

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog_v2[]> Logs { get { return logs; } }

        public IChildAnalyzer ChildAnalyzer { get; set; }

        public void BuildMatrices()
        {
            if (ChildAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");
            if (currentSimulation < 0)
                InitializeCoefficientsProvider();

            provider.Reset();
            coefficientsProvider.CurrentOrder = -1;
            ChildAnalyzer.BuildMatrices();

            matrices[0] = new Dictionary<int, IMatrix>(linearSystems.Count);
            foreach (var linearSystem in linearSystems.Values)
            {
                matrices[0].Add(linearSystem.Subdomain.ID, linearSystem.Matrix.Copy(false));
            }
            for (int i = 0; i < expansionOrder; i++)
            {
                provider.Reset();
                coefficientsProvider.CurrentOrder = i;
                ChildAnalyzer.BuildMatrices();

                matrices[i + 1] = new Dictionary<int, IMatrix>(linearSystems.Count);
                foreach (var linearSystem in linearSystems.Values)
                {
                    matrices[i + 1].Add(linearSystem.Subdomain.ID, linearSystem.Matrix.Copy(false));
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
            var matricesPerSubdomain = new Dictionary<int, IMatrix[]>();
            foreach (var linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                var tempMatrices = new IMatrix[expansionOrder + 1];
                for (int i = 0; i <= expansionOrder; i++)
                    tempMatrices[i] = matrices[i][id];
                matricesPerSubdomain.Add(id, tempMatrices);
            }

            foreach (var linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                IMatrix combo = matricesPerSubdomain[id][0].Scale(coefficients[0]);
                for (int i = 1; i < coefficients.Count; ++i) combo.AxpyIntoThis(matricesPerSubdomain[id][i], coefficients[i]);
                linearSystem.Matrix = combo;
            }

            //for (int i = 0; i < expansionOrder; i++)
            //{
            //    foreach (var subdomain in subdomains.Values)
            //    {
            //        SkylineMatrix2D<double> k = (SkylineMatrix2D<double>)matrices[i + 1][subdomain.ID];
            //        subdomain.Matrix.LinearCombination(new double[] { randomNumbers[currentSimulation][i] }, new IMatrix2D<double>[] { k });
            //    }
            //}
        }

        public void Initialize(bool isFirstAnalysis)
        {
            if (ChildAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");
            ChildAnalyzer.Initialize(isFirstAnalysis);
        }

        public void Solve()
        {
            if (ChildAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");

            string[] values = new string[simulations];
            for (int i = 0; i < simulations; i++)
            {
                currentSimulation = i;
                ComposeStochasticMatrixFromMatrices();
                ChildAnalyzer.Initialize(false);
                ChildAnalyzer.Solve();
                values[i] = linearSystems[1].Solution[28].ToString();
            }

            File.WriteAllLines(String.Format(@"montecarlo{0}.txt", expansionOrder), values);
        }

        #endregion
    }
}
