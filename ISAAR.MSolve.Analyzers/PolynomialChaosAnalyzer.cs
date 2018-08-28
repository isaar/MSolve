using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using Troschuetz.Random.Distributions.Continuous;
using System.IO;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.Analyzers
{
    public class PolynomialChaosAnalyzer : IAnalyzer
    {
        private int currentSimulation = -1;
        private readonly bool shouldFactorizeMatrices;
        private readonly int expansionOrder;
        private readonly int simulations;
        private readonly IDictionary<int, ILinearSystem> subdomains;
        private readonly IDictionary<int, IMatrix2D>[] matrices;
        private readonly IList<Dictionary<int, SkylineMatrix2D>> factorizedMatrices = new List<Dictionary<int, SkylineMatrix2D>>();
        private readonly Model model;
        private readonly Dictionary<int, IAnalyzerLog[]> logs = new Dictionary<int, IAnalyzerLog[]>();
        private readonly IAnalyzerProvider provider;
        private readonly double[][] randomNumbers;
        private IAnalyzer childAnalyzer;
        private IAnalyzer parentAnalyzer = null;
        //private readonly LognormalPCFileStochasticCoefficientsProvider coefficientsProvider;
        private readonly IPCCoefficientsProvider coefficientsProvider;
        private List<int>[] nonZeroPsi;

        public IDictionary<int, IMatrix2D>[] Matrices { get { return matrices; } }
        public IList<Dictionary<int, SkylineMatrix2D>> FactorizedMatrices { get { return factorizedMatrices; } }

        public PolynomialChaosAnalyzer(Model model, IAnalyzerProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ILinearSystem> subdomains, IPCCoefficientsProvider coefficientsProvider, int expansionOrder, int simulations, bool shouldFactorizeMatrices)
        {
            this.shouldFactorizeMatrices = shouldFactorizeMatrices;
            this.childAnalyzer = embeddedAnalyzer;
            this.provider = provider;
            this.model = model;
            this.subdomains = subdomains;
            this.expansionOrder = coefficientsProvider.ExpansionOrder;
            this.simulations = simulations;
            this.childAnalyzer.ParentAnalyzer = this;
            this.matrices = new Dictionary<int, IMatrix2D>[coefficientsProvider.NoOfMatrices + 1];
            this.randomNumbers = new double[simulations][];
            this.coefficientsProvider = coefficientsProvider;

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

        public PolynomialChaosAnalyzer(Model model, IAnalyzerProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ILinearSystem> subdomains, IPCCoefficientsProvider coefficientsProvider, int expansionOrder, int simulations)
            : this(model, provider, embeddedAnalyzer, subdomains, coefficientsProvider, expansionOrder, simulations, true)
        {
        }

        private void InitializeCoefficientsProvider()
        {
            //coefficientsProvider = new FileStochasticCoefficientsProvider("Lognormal.csv", stochasticDomain);
            //coefficientsProvider = new FileStochasticCoefficientsProvider("Gaussian.txt", 50, '\t', stochasticDomain);
            foreach (var subdomain in model.Subdomains)
                foreach (var e in subdomain.ElementsDictionary.Values.Where(e => e.ElementType is IStochasticFiniteElement))
                    ((IStochasticFiniteElement)e.ElementType).CoefficientsProvider = coefficientsProvider;
        }

        private void MakePreconditioners()
        {
            factorizedMatrices.Clear();
            for (int i = 0; i < matrices.Length; i++)
                factorizedMatrices.Add(new Dictionary<int, SkylineMatrix2D>());

            foreach (var sub in subdomains)
            {
                //for (int i = 0; i < matrices.Length; i++)
                for (int i = 0; i < 1; i++)
                {
                    var m = (SkylineMatrix2D)((SkylineMatrix2D)matrices[i][sub.Key]).Clone();
                    m.Factorize(1e-32, new List<IVector>(), new List<int>());

                    if (factorizedMatrices[i].ContainsKey(sub.Key))
                        factorizedMatrices[i][sub.Key] = m;
                    else
                        factorizedMatrices[i].Add(sub.Key, m);
                    //if (factorizedMatrices.Count == 0)
                    //    factorizedMatrices.Add(new Dictionary<int, SkylineMatrix2D<double>>());
                    //if (factorizedMatrices[0].ContainsKey(sub.Key))
                    //    factorizedMatrices[0][sub.Key] = m;
                    //else
                    //    factorizedMatrices[0].Add(sub.Key, m);
                }
            }
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
            if (childAnalyzer == null) throw new InvalidOperationException("Polynomial chaos analyzer must contain an embedded analyzer.");
            if (currentSimulation < 0)
                InitializeCoefficientsProvider();

            for (int i = 0; i < coefficientsProvider.NoOfMatrices; i++)
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

            provider.Reset();
            coefficientsProvider.CurrentOrder = -1;
            childAnalyzer.BuildMatrices();

            matrices[0] = new Dictionary<int, IMatrix2D>(subdomains.Count);
            foreach (var subdomain in subdomains.Values)
            {
                SkylineMatrix2D k = (SkylineMatrix2D)subdomain.Matrix;
                matrices[0].Add(subdomain.ID, (SkylineMatrix2D)k.Clone());
            }
        }

        //private void BuildStochasticMatrices()
        //{
        //    var currentRandomNumbers = randomNumbers[currentSimulation];
        //    var coefficients = new double[] { 1 }.Concat(currentRandomNumbers).ToList<double>();
        //    var matricesPerSubdomain = new Dictionary<int, IMatrix2D<double>[]>();
        //    foreach (var subdomain in subdomains.Values)
        //    {
        //        int id = subdomain.ID;
        //        var tempMatrices = new IMatrix2D<double>[expansionOrder + 1];
        //        for (int i = 0; i <= expansionOrder; i++)
        //            tempMatrices[i] = matrices[i][id];
        //        matricesPerSubdomain.Add(id, tempMatrices);
        //    }

        //    //foreach (var subdomain in subdomains.Values)
        //    //    subdomain.Matrix = (SkylineMatrix2D<double>)((SkylineMatrix2D<double>)matrices[0][subdomain.ID]).Clone();
        //    //foreach (var subdomain in subdomains.Values)
        //    //    subdomain.Matrix.LinearCombination(coefficients, matricesPerSubdomain[subdomain.ID]);
        //}

        public void Initialize()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Polynomial chaos analyzer must contain an embedded analyzer.");
        }

        private void InitializeNonZeroPsi()
        {
            const double tolerance = 1e-15;
            int psiSize = coefficientsProvider.Calculator.PsiSize;
            nonZeroPsi = new List<int>[psiSize];
            for (int i = 0; i < psiSize; i++)
            {
                var nonZeros = new List<int>();
                var currentPsiBasis = coefficientsProvider.Calculator.PsiBasis[i];
                for (int j = 0; j < currentPsiBasis.Length; j++)
                    if (Math.Abs(currentPsiBasis[j]) > tolerance)
                        nonZeros.Add(j);

                nonZeroPsi[i] = nonZeros;
            }
        }

        public void Solve()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Polynomial chaos analyzer must contain an embedded analyzer.");

            //BuildStochasticMatrices();
            InitializeNonZeroPsi();
            if (shouldFactorizeMatrices) MakePreconditioners();
            childAnalyzer.Initialize();
            childAnalyzer.Solve();

            int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[1][DOFType.X];
            //int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[84][DOFType.Y];
            //int dofNo = 112;
            var subdomain = subdomains.Select(x => x.Value).First();
            int meanLength = subdomain.Solution.Length / coefficientsProvider.Calculator.PsiSize;
            double[] values = new double[simulations];
            string[] valuesString = new string[simulations];
            for (int m = 0; m < simulations; m++)
            {
                var psii = EvaluatePSI(randomNumbers[m]);
                for (int i = 0; i < coefficientsProvider.Calculator.PsiSize; i++)
                    values[m] += subdomain.Solution[i * meanLength + dofNo] * psii[i];
                valuesString[m] = values[m].ToString();
            }
            File.WriteAllLines(String.Format(@"polychaos{0}-{1}-{2}.txt", coefficientsProvider.OrderM, coefficientsProvider.OrderP, expansionOrder), valuesString);

            //string[] values = new string[simulations];
            //for (int i = 0; i < simulations; i++)
            //{
            //    currentSimulation = i;
            //    values[i] = subdomains[1].Solution[28].ToString();
            //}

            //File.WriteAllLines(String.Format(@"polychaos{0}.txt", expansionOrder), values);
        }

        #endregion

        private double EvaluateHermite(int p, double ksi)
        {
            double result = 0;
            var hp = coefficientsProvider.Calculator.HermitePolynomials[p - 1];
            for (int k = 0; k < hp.Length; k++)
                result += hp[k] * Math.Pow(ksi, p - k);

            return result;
        }
        //function[Her] = EvaluateHermite(Hermitefunctions,p,ksi)
        //for i = p
        //    Her = 0;
        //        for k = 1:size(Hermitefunctions{i},2)
        //            t = Hermitefunctions{i}(k);
        //            Her = Her + t*ksi^(i-(k-1)) ;
        //        end
        //    end
        //end

        private double[] EvaluatePSI(double[] ksi)
        {
            var psiSize = coefficientsProvider.Calculator.PsiSize;
            var psii = new double[psiSize];
            for (int i = 0; i < psiSize; i++)
            {
                psii[i] = 1;
                for (int j = 0; j < nonZeroPsi[i].Count; j++)
                {
                    var index = nonZeroPsi[i][j];
                    var x = ksi[index];
                    var p = coefficientsProvider.Calculator.PsiBasis[i][index];
                    var hermite = EvaluateHermite(p, x);
                    psii[i] *= hermite;
                }
            }

            return psii;
        }

        //function [Psi] = EvaluatePsi(PC,ksi);
        // PsiSize = size(PC.PC.PsiBasis,1);
        //for i = 1:PsiSize
        //    Psi{i,1} = 1;
        //    ind = find(PC.PC.PsiBasis{i});
        //    if size(ind,2) == 0
        //    else
        //    for j = 1:size(ind,2)
        //        x = ksi(ind(j));
        //        p =  PC.PC.PsiBasis{i}(ind(j));
        //        Her = EvaluateHermite(PC.Hermite,p,x);
        //        Psi{i} = Psi{i,1}*Her;
        //    end
        //end
    }
}
