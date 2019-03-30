using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Input;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;
using Troschuetz.Random.Distributions.Continuous;

namespace ISAAR.MSolve.Analyzers
{
    public enum StiffnessMatrixProductionMode
    {
        Normal = 0,
        StoreToDisk,
        LoadFromDiskAndCalculate,
        StoreToDiskAndCalculate
    }

    public class MonteCarloAnalyzerWithStochasticMaterial_v2 : IParentAnalyzer
    {
        private int currentSimulation = -1;
        private readonly int blockSize = 5;
        private readonly int expansionOrder;
        private readonly int simulations;
        private readonly int simulationStartFrom = 0;
        private readonly int randomFileSimulations = 50000;
        private readonly IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems;
        //private readonly IDictionary<int, IMatrix<double>> matrices;
        private readonly IDictionary<int, IMatrix>[] matrices;
        private readonly Model_v2 model;
        private readonly IAnalyzerProvider_v2 provider;
        private double[][] randomNumbers;
        private readonly StiffnessMatrixProductionMode stiffnessMatrixProductionMode = StiffnessMatrixProductionMode.Normal;
        private readonly string stiffnessMatrixPath = String.Empty;
        private readonly string randomsReadFileName = String.Empty;
        private readonly string fileNameForLogging = "monteCarlo";
        private readonly IStochasticMaterialCoefficientsProvider coefficientsProvider;
        private readonly List<int> matrixOrder = new List<int>();
        private readonly List<double> matrixMagnitudes = new List<double>();

        public IDictionary<int, LdlSkyline> FactorizedMatrices { get; } = new Dictionary<int, LdlSkyline>();

        public MonteCarloAnalyzerWithStochasticMaterial_v2(Model_v2 model, IAnalyzerProvider_v2 provider, 
            IChildAnalyzer embeddedAnalyzer, ISolver_v2 solver, IStochasticMaterialCoefficientsProvider coefficientsProvider, 
            int expansionOrder, int simulations)
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
            this.coefficientsProvider = coefficientsProvider;
            //this.stochasticDomain = stochasticDomain;
        }

        public MonteCarloAnalyzerWithStochasticMaterial_v2(Model_v2 model, IAnalyzerProvider_v2 provider, 
            IChildAnalyzer embeddedAnalyzer, ISolver_v2 solver, IStochasticMaterialCoefficientsProvider coefficientsProvider,
            int expansionOrder, int simulations, string fileNameForLogging)
            : this(model, provider, embeddedAnalyzer, solver, coefficientsProvider, expansionOrder, simulations)
        {
            this.fileNameForLogging = fileNameForLogging;
        }

        public MonteCarloAnalyzerWithStochasticMaterial_v2(Model_v2 model, IAnalyzerProvider_v2 provider, 
            IChildAnalyzer embeddedAnalyzer, ISolver_v2 solver, IStochasticMaterialCoefficientsProvider coefficientsProvider,
            int expansionOrder, int simulations, StiffnessMatrixProductionMode stiffnessMatrixProductionMode, 
            string fileNameForLogging, string stiffnessMatrixPath)
            : this(model, provider, embeddedAnalyzer, solver, coefficientsProvider, expansionOrder, simulations, 
                  fileNameForLogging)
        {
            this.stiffnessMatrixPath = stiffnessMatrixPath;
            this.stiffnessMatrixProductionMode = stiffnessMatrixProductionMode;
        }

        public MonteCarloAnalyzerWithStochasticMaterial_v2(Model_v2 model, IAnalyzerProvider_v2 provider, 
            IChildAnalyzer embeddedAnalyzer, ISolver_v2 solver, IStochasticMaterialCoefficientsProvider coefficientsProvider,
            int expansionOrder, int simulations, int blockSize, StiffnessMatrixProductionMode stiffnessMatrixProductionMode, 
            string fileNameForLogging, string stiffnessMatrixPath)
            : this(model, provider, embeddedAnalyzer, solver, coefficientsProvider, expansionOrder, simulations, 
                  stiffnessMatrixProductionMode, fileNameForLogging, stiffnessMatrixPath)
        {
            this.blockSize = blockSize;
        }

        public MonteCarloAnalyzerWithStochasticMaterial_v2(Model_v2 model, IAnalyzerProvider_v2 provider, 
            IChildAnalyzer embeddedAnalyzer, ISolver_v2 solver, IStochasticMaterialCoefficientsProvider coefficientsProvider,
            int expansionOrder, int simulations, int blockSize, StiffnessMatrixProductionMode stiffnessMatrixProductionMode, 
            string fileNameForLogging, string stiffnessMatrixPath, string randomsReadFileName, 
            int simulationStartFrom)
            : this(model, provider, embeddedAnalyzer, solver, coefficientsProvider, expansionOrder, simulations, blockSize, stiffnessMatrixProductionMode, fileNameForLogging, stiffnessMatrixPath)
        {
            this.randomsReadFileName = randomsReadFileName;
            this.simulationStartFrom = simulationStartFrom;
        }

        //private void InitializeCoefficientsProvider()
        //{
        //    //coefficientsProvider = new FileStochasticCoefficientsProvider("Lognormal.csv", stochasticDomain);
        //    //coefficientsProvider = new FileStochasticCoefficientsProvider("Gaussian.txt", 50, '\t', stochasticDomain);
        //    foreach (var subdomain in model.Subdomains)
        //        foreach (var e in subdomain.ElementsDictionary.Values.Where(e => e.ElementType is IStochasticFiniteElement))
        //            ((IStochasticFiniteElement)e.ElementType).CoefficientsProvider = coefficientsProvider;
        //}

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog_v2[]> Logs { get; } = new Dictionary<int, IAnalyzerLog_v2[]>();

        public IChildAnalyzer ChildAnalyzer { get; set; }

        public void BuildMatrices()
        {
            if (ChildAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");
            if (currentSimulation < 0)
            {
                if (stiffnessMatrixProductionMode == StiffnessMatrixProductionMode.LoadFromDiskAndCalculate)
                {
                    var fileName = String.Format(@"{0}-Order.txt", fileNameForLogging);
                    using (var sw = File.OpenText(fileName))
                    {
                        while (!sw.EndOfStream)
                            matrixOrder.Add(Int32.Parse(sw.ReadLine()));
                    }
                }
                return;
            }

            if (stiffnessMatrixProductionMode == StiffnessMatrixProductionMode.LoadFromDiskAndCalculate)
                ReadMatricesFromFile(currentSimulation);
            else
            {
                provider.Reset();
                ChildAnalyzer.BuildMatrices();
            }
            WriteMatricesToFile(currentSimulation);
        }

        //private void ComposeStochasticMatrixFromMatrices()
        //{
        //    var currentRandomNumbers = randomNumbers[currentSimulation];
        //    var coefficients = new double[] { 1 }.Concat(currentRandomNumbers).ToList<double>();
        //    var matricesPerSubdomain = new Dictionary<int, IMatrix<double>[]>();
        //    foreach (var subdomain in subdomains.Values)
        //    {
        //        int id = subdomain.ID;
        //        var tempMatrices = new IMatrix<double>[expansionOrder + 1];
        //        for (int i = 0; i <= expansionOrder; i++)
        //            tempMatrices[i] = matrices[i][id];
        //        matricesPerSubdomain.Add(id, tempMatrices);
        //    }

        //    foreach (var subdomain in subdomains.Values)
        //        subdomain.Matrix = (SkylineMatrix<double>)((SkylineMatrix<double>)matrices[0][subdomain.ID]).Clone();
        //    foreach (var subdomain in subdomains.Values)
        //        subdomain.Matrix.LinearCombination(coefficients, matricesPerSubdomain[subdomain.ID]);

            
        //    //for (int i = 0; i < expansionOrder; i++)
        //    //{
        //    //    foreach (var subdomain in subdomains.Values)
        //    //    {
        //    //        SkylineMatrix<double> k = (SkylineMatrix<double>)matrices[i + 1][subdomain.ID];
        //    //        subdomain.Matrix.LinearCombination(new double[] { randomNumbers[currentSimulation][i] }, new IMatrix<double>[] { k });
        //    //    }
        //    //}
        //}

        public void Initialize(bool isFirstAnalysis = true)
        {
            if (ChildAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");

            if (String.IsNullOrEmpty(randomsReadFileName))
            {
                randomNumbers = new double[simulations][];
                NormalDistribution n = new NormalDistribution();
                n.Mu = 0;
                n.Sigma = 1;
                string[] randoms = new string[simulations];
                for (int i = 0; i < simulations; i++)
                {
                    randomNumbers[i] = new double[expansionOrder];
                    for (int j = 0; j < expansionOrder; j++)
                        randomNumbers[i][j] = n.NextDouble();
                }
                //using (var sw = File.CreateText(String.Format(@"randoms{0}.txt", expansionOrder)))
                //{
                //    for (int j = 0; j < expansionOrder; j++)
                //        for (int i = 0; i < simulations; i++)
                //            sw.WriteLine(randomNumbers[i][j]);
                //}
            }
            else
            {
                randomNumbers = new double[randomFileSimulations][];
                for (int i = 0; i < randomFileSimulations; i++)
                    randomNumbers[i] = new double[expansionOrder];
                using (var sw = File.OpenText(randomsReadFileName))
                {
                    for (int j = 0; j < expansionOrder; j++)
                        for (int i = 0; i < randomFileSimulations; i++)
                            randomNumbers[i][j] = Double.Parse(sw.ReadLine());
                }
            }

            if (stiffnessMatrixProductionMode != StiffnessMatrixProductionMode.LoadFromDiskAndCalculate) return;
        }

        public void Solve()
        {
            if (ChildAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");

            if (stiffnessMatrixProductionMode == StiffnessMatrixProductionMode.LoadFromDiskAndCalculate)
                SolveWithOrder();
            else
                SolveNormal();
        }

        private void SolveNormal()
        {
            //int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[150][DOFType.Y];
            //int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[84][DOFType.Y];
            int dofNo = model.Subdomains[0].FreeDofOrdering.FreeDofs[model.NodesDictionary[10], DOFType.Y];
            //int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[450][DOFType.Y];
            //int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[601][DOFType.Y];
            //int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[6051][DOFType.Y];
            //int dofNo = model.Subdomains[36].NodalDOFsDictionary[6051][DOFType.Y];
            string[] values = new string[simulations];
            double[] numberValues = new double[simulations];
            var fileName = String.Format(@"{0}-{1}-{2}.txt", fileNameForLogging, expansionOrder, simulationStartFrom);
            var fileNameIterations = String.Format(@"{0}-{1}-{2}-Iters.txt", fileNameForLogging, expansionOrder, simulationStartFrom);
            var fileNameTimes = String.Format(@"{0}-{1}-{2}-Times.txt", fileNameForLogging, expansionOrder, simulationStartFrom);
            StreamWriter sw;
            sw = File.CreateText(fileName);
            sw.Dispose();
            sw = File.CreateText(fileNameIterations);
            sw.Dispose();
            var times = new Dictionary<string, TimeSpan>();
            times.Add("all", TimeSpan.Zero);
            times.Add("element", TimeSpan.Zero);
            times.Add("factorize", TimeSpan.Zero);
            times.Add("solution", TimeSpan.Zero);
            var start = DateTime.Now;
            var iterationCount = new int[simulations - simulationStartFrom];
            int totalIterations = 0;
            for (int i = simulationStartFrom; i < simulationStartFrom + simulations; i++)
            {
                currentSimulation = i;
                coefficientsProvider.RandomVariables = randomNumbers[currentSimulation];
                var e = DateTime.Now;
                BuildMatrices();
                times["element"] += DateTime.Now - e;
                if (stiffnessMatrixProductionMode == StiffnessMatrixProductionMode.StoreToDisk) continue;

                e = DateTime.Now;
                ChildAnalyzer.Initialize(false);
                times["factorize"] += DateTime.Now - e;
                GCSettings.LatencyMode = GCLatencyMode.LowLatency;
                e = DateTime.Now;
                ChildAnalyzer.Solve();

                times["solution"] += DateTime.Now - e;
                GCSettings.LatencyMode = GCLatencyMode.Batch;
                values[i - simulationStartFrom] = linearSystems[0].Solution[dofNo].ToString();
                numberValues[i - simulationStartFrom] = linearSystems[0].Solution[dofNo];
                //values[i] = subdomains[36].Solution[dofNo].ToString();

                //values[i] = matrixMagnitudes[i].ToString();
                using (sw = File.AppendText(fileName))
                    sw.WriteLine(values[i - simulationStartFrom]);
                using (sw = File.AppendText(fileNameIterations))
                    sw.WriteLine(iterationCount[i - simulationStartFrom].ToString());
                //using (sw = File.CreateText(fileNameTimes))
                //{
                //    sw.WriteLine(String.Format("Elements: {0}", times["element"].ToString()));
                //    sw.WriteLine(String.Format("Factorize: {0}", times["factorize"].ToString()));
                //    sw.WriteLine(String.Format("Solution: {0}", times["solution"].ToString()));
                //    sw.WriteLine(String.Format("All: {0}", times["all"].ToString()));
                //    sw.WriteLine(String.Format("Total iterations: {0}", totalIterations.ToString()));
                //}
            }
            MonteCarloMeanValue = numberValues.Average();
            double sumOfSquaresOfDifferences = numberValues.Select(val => (val - MonteCarloMeanValue) * (val - MonteCarloMeanValue)).Sum();
            MonteCarloStandardDeviation = Math.Sqrt(sumOfSquaresOfDifferences / numberValues.Length);
            times["all"] = DateTime.Now - start;

            //File.WriteAllLines(String.Format(@"{0}-{1}.txt", fileNameForLogging, expansionOrder), values);
        }

        public double MonteCarloMeanValue { get; set; }
        public double MonteCarloStandardDeviation { get; set; }

        private void MakePreconditioner(int simulation)
        {
            int matrixNo = matrixOrder[simulation + blockSize / 2];
            string name = String.IsNullOrWhiteSpace(stiffnessMatrixPath) ? "K" : stiffnessMatrixPath;
            string path = Path.GetDirectoryName(name);
            string nameOnly = Path.GetFileNameWithoutExtension(name);
            string ext = Path.GetExtension(name);

            foreach (var linearSystem in linearSystems)
            {
                SkylineMatrix m = SkylineMatrixReader.ReadFromSimilarlyNamedFiles(
                    String.Format("{0}\\{1}Sub{3}Sim{4}{2}", path, nameOnly, ext, linearSystem.Key, matrixNo));
                LdlSkyline factor = m.FactorLdl(true, 1e-8);
                if (FactorizedMatrices.ContainsKey(linearSystem.Key)) FactorizedMatrices[linearSystem.Key] = factor;
                else FactorizedMatrices.Add(linearSystem.Key, factor);
            }
        }

        private void SolveWithOrder()
        {
            int dofNo = model.Subdomains[0].FreeDofOrdering.FreeDofs[model.NodesDictionary[6051], DOFType.Y];
            string[] values = new string[simulations];
            var fileName = String.Format(@"{0}-{1}-{2}.txt", fileNameForLogging, expansionOrder, simulationStartFrom);
            StreamWriter sw = File.CreateText(fileName);
            sw.Dispose();
            var times = new Dictionary<string, TimeSpan>();
            times.Add("all", TimeSpan.Zero);
            times.Add("element", TimeSpan.Zero);
            times.Add("factorize", TimeSpan.Zero);
            times.Add("solution", TimeSpan.Zero);
            var start = DateTime.Now;
            var e = start;
            for (int i = simulationStartFrom; i < simulationStartFrom + simulations; i++)
            {
                if (i % blockSize == 0)
                {
                    e = DateTime.Now;
                    MakePreconditioner(i);
                    times["factorize"] += DateTime.Now - e;
                }
                currentSimulation = matrixOrder[i];
                //coefficientsProvider.RandomVariables = randomNumbers[currentSimulation];
                e = DateTime.Now;
                BuildMatrices();
                times["element"] += DateTime.Now - e;
                if (stiffnessMatrixProductionMode == StiffnessMatrixProductionMode.StoreToDisk) continue;

                e = DateTime.Now;
                ChildAnalyzer.Initialize(false);
                ChildAnalyzer.Solve();
                times["solution"] += DateTime.Now - e;
                values[i] = linearSystems[1].Solution[dofNo].ToString();

                //using (sw = File.AppendText(fileName))
                //{
                //    sw.WriteLine(values[i]);
                //}

            }
            times["all"] = DateTime.Now - start;
            var s = new List<string>();
            s.Add(times["all"].ToString());
            s.Add(times["element"].ToString());
            s.Add(times["factorize"].ToString());
            s.Add(times["solution"].ToString());
            File.WriteAllLines(String.Format(@"{0}-Times-{1}-{2}.txt", fileNameForLogging, blockSize, simulationStartFrom), s);
        }

        private void WriteMatricesToFile(int simulation)
        {
            if (stiffnessMatrixProductionMode == StiffnessMatrixProductionMode.Normal || stiffnessMatrixProductionMode == StiffnessMatrixProductionMode.LoadFromDiskAndCalculate) return;
            string name = String.IsNullOrWhiteSpace(stiffnessMatrixPath) ? "K" : stiffnessMatrixPath;

            string path = Path.GetDirectoryName(name);
            string nameOnly = Path.GetFileNameWithoutExtension(name);
            string ext = Path.GetExtension(name);

            //foreach (var sub in subdomains)
            //    sub.Value.Matrix.WriteToFile(String.Format(@"{0}\{1}Sub{3}Sim{4}{2}", path, nameOnly, ext, sub.Key, simulation));
        }

        private void ReadMatricesFromFile(int simulation)
        {
            if (stiffnessMatrixProductionMode != StiffnessMatrixProductionMode.LoadFromDiskAndCalculate) return;
            string name = String.IsNullOrWhiteSpace(stiffnessMatrixPath) ? "K" : stiffnessMatrixPath;

            string path = Path.GetDirectoryName(name);
            string nameOnly = Path.GetFileNameWithoutExtension(name);
            string ext = Path.GetExtension(name);

            foreach (var linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                SkylineMatrix m = SkylineMatrixReader.ReadFromSimilarlyNamedFiles(
                    String.Format("{0}\\{1}Sub{3}Sim{4}{2}", path, nameOnly, ext, id, simulation));
                linearSystem.Matrix = m;

                //double d = 0;
                //for (int i = 0; i < m.RowIndex.Length - 1; i++)
                //    d += m.Data[m.RowIndex[i]];
                //matrixMagnitudes.Add(d / (double)m.RowIndex.Length);
            }
        }

        #endregion
    }
}
