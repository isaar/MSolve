using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using Troschuetz.Random.Distributions.Continuous;
using System.IO;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    public class VRFMonteCarloAnalyzerWithStochasticMaterial : IAnalyzer
    {
        private int currentSimulation = -1;
        private readonly int blockSize = 5;
        private readonly int expansionOrder;
        private readonly int simulations;
        private readonly int simulationStartFrom = 0;
        private readonly int randomFileSimulations = 50000;
        private readonly IDictionary<int, ILinearSystem> subdomains;
        //private readonly IDictionary<int, IMatrix2D<double>> matrices;
        private readonly IDictionary<int, IMatrix2D>[] matrices;
        private readonly IDictionary<int, SkylineMatrix2D> factorizedMatrices = new Dictionary<int, SkylineMatrix2D>();
        private readonly Model model;
        private readonly Dictionary<int, IAnalyzerLog[]> logs = new Dictionary<int, IAnalyzerLog[]>();
        private readonly IAnalyzerProvider provider;
        private double[][] randomNumbers;
        //private readonly StiffnessMatrixProductionMode stiffnessMatrixProductionMode = StiffnessMatrixProductionMode.Normal;
        private readonly string stiffnessMatrixPath = String.Empty;
        private readonly string randomsReadFileName = String.Empty;
        //private readonly double[] stochasticDomain;
        private IAnalyzer childAnalyzer;
        private IAnalyzer parentAnalyzer = null;
        private readonly string fileNameForLogging = "monteCarlo";
        private readonly PowerSpectrumTargetEvaluatorCoefficientsProvider coefficientsProvider;
        private readonly IStochasticMaterialCoefficientsProvider approximateCoefficientsProvider; 
        private readonly List<int> matrixOrder = new List<int>();
        private readonly List<double> matrixMagnitudes = new List<double>();
        private Dictionary<int, double[]> vrfs;

        public IDictionary<int, SkylineMatrix2D> FactorizedMatrices { get { return factorizedMatrices; } }

        public VRFMonteCarloAnalyzerWithStochasticMaterial(Model model, IAnalyzerProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ILinearSystem> subdomains, PowerSpectrumTargetEvaluatorCoefficientsProvider coefficientsProvider, IStochasticMaterialCoefficientsProvider approximateCoefficientsProvider, int expansionOrder, int simulations)
        {
            this.childAnalyzer = embeddedAnalyzer;
            this.provider = provider;
            this.model = model;
            this.subdomains = subdomains;
            this.expansionOrder = expansionOrder;
            this.childAnalyzer.ParentAnalyzer = this;
            //this.matrices = new Dictionary<int, IMatrix2D<double>>(subdomains.Count);
            this.matrices = new Dictionary<int, IMatrix2D>[expansionOrder + 1];
            this.coefficientsProvider = coefficientsProvider;
            this.approximateCoefficientsProvider = approximateCoefficientsProvider;
            this.simulations = coefficientsProvider.NPhi * coefficientsProvider.NPtsVRF;
            //this.stochasticDomain = stochasticDomain;
        }

        public VRFMonteCarloAnalyzerWithStochasticMaterial(Model model, IAnalyzerProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ILinearSystem> subdomains, PowerSpectrumTargetEvaluatorCoefficientsProvider coefficientsProvider, IStochasticMaterialCoefficientsProvider approximateCoefficientsProvider,
            int expansionOrder, int simulations, string fileNameForLogging)
            : this(model, provider, embeddedAnalyzer, subdomains, coefficientsProvider, approximateCoefficientsProvider, expansionOrder, simulations)
        {
            this.fileNameForLogging = fileNameForLogging;
        }

        public VRFMonteCarloAnalyzerWithStochasticMaterial(Model model, IAnalyzerProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ILinearSystem> subdomains, PowerSpectrumTargetEvaluatorCoefficientsProvider coefficientsProvider, IStochasticMaterialCoefficientsProvider approximateCoefficientsProvider,
            int expansionOrder, int simulations, StiffnessMatrixProductionMode stiffnessMatrixProductionMode, string fileNameForLogging, string stiffnessMatrixPath)
            : this(model, provider, embeddedAnalyzer, subdomains, coefficientsProvider, approximateCoefficientsProvider, expansionOrder, simulations, fileNameForLogging)
        {
            this.stiffnessMatrixPath = stiffnessMatrixPath;
            //this.stiffnessMatrixProductionMode = stiffnessMatrixProductionMode;
        }

        public VRFMonteCarloAnalyzerWithStochasticMaterial(Model model, IAnalyzerProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ILinearSystem> subdomains, PowerSpectrumTargetEvaluatorCoefficientsProvider coefficientsProvider, IStochasticMaterialCoefficientsProvider approximateCoefficientsProvider,
            int expansionOrder, int simulations, int blockSize, StiffnessMatrixProductionMode stiffnessMatrixProductionMode, string fileNameForLogging, string stiffnessMatrixPath)
            : this(model, provider, embeddedAnalyzer, subdomains, coefficientsProvider, approximateCoefficientsProvider, expansionOrder, simulations, stiffnessMatrixProductionMode, fileNameForLogging, stiffnessMatrixPath)
        {
            this.blockSize = blockSize;
        }

        public VRFMonteCarloAnalyzerWithStochasticMaterial(Model model, IAnalyzerProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ILinearSystem> subdomains, PowerSpectrumTargetEvaluatorCoefficientsProvider coefficientsProvider, IStochasticMaterialCoefficientsProvider approximateCoefficientsProvider,
            int expansionOrder, int simulations, int blockSize, StiffnessMatrixProductionMode stiffnessMatrixProductionMode, string fileNameForLogging, string stiffnessMatrixPath, string randomsReadFileName, 
            int simulationStartFrom)
            : this(model, provider, embeddedAnalyzer, subdomains, coefficientsProvider, approximateCoefficientsProvider, expansionOrder, simulations, blockSize, stiffnessMatrixProductionMode, fileNameForLogging, stiffnessMatrixPath)
        {
            this.randomsReadFileName = randomsReadFileName;
            this.simulationStartFrom = simulationStartFrom;
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
            provider.Reset();
            childAnalyzer.BuildMatrices();
        }

        public void Initialize()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");

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
        }

        public void Solve()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");

            SolveNormal();
        }

        private void SolveNormal()
        {
            //int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[150][DOFType.Y];
            //int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[84][DOFType.Y];
            int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[80][DOFType.X];
            //int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[450][DOFType.Y];
            //int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[601][DOFType.Y];
            //int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[6051][DOFType.Y];
            //int dofNo = model.Subdomains[36].NodalDOFsDictionary[6051][DOFType.Y];
            string[] values = new string[simulations];
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
            currentSimulation = 0;
            var dofValuesPerPhi = new Dictionary<int, double[][]>(model.Nodes.Count);
            foreach (var d in model.NodalDOFsDictionary)
                dofValuesPerPhi.Add(d.Key, new double[][] { new double[coefficientsProvider.NPhi], new double[coefficientsProvider.NPhi], new double[coefficientsProvider.NPhi] });
            var vrfPerNPtsVRF = new Dictionary<int, double[][]>(model.Nodes.Count);
            foreach (var d in model.NodalDOFsDictionary)
                vrfPerNPtsVRF.Add(d.Key, new double[][] { new double[coefficientsProvider.NPtsVRF], new double[coefficientsProvider.NPtsVRF], new double[coefficientsProvider.NPtsVRF] });

            for (int i = 0; i < coefficientsProvider.NPtsVRF; i++)
            {
                coefficientsProvider.CurrentFrequency = i;
                for (int j = 0; j < coefficientsProvider.NPhi; j++)
                {
                    coefficientsProvider.RandomVariables = randomNumbers[currentSimulation];
                    coefficientsProvider.CurrentMCS = j;
                    var e = DateTime.Now;
                    BuildMatrices();
                    times["element"] += DateTime.Now - e;
                    //if (stiffnessMatrixProductionMode == StiffnessMatrixProductionMode.StoreToDisk) continue;

                    e = DateTime.Now;
                    childAnalyzer.Initialize();
                    times["factorize"] += DateTime.Now - e;
                    GCSettings.LatencyMode = GCLatencyMode.LowLatency;
                    e = DateTime.Now;
                    childAnalyzer.Solve();

                    times["solution"] += DateTime.Now - e;
                    GCSettings.LatencyMode = GCLatencyMode.Batch;
                    values[i - simulationStartFrom] = subdomains[1].Solution[dofNo].ToString();
                    //values[i] = subdomains[36].Solution[dofNo].ToString();

                    //values[i] = matrixMagnitudes[i].ToString();
                    using (sw = File.AppendText(fileName))
                        //sw.WriteLine(values[i - simulationStartFrom]);
                        sw.WriteLine(GetCSVText(subdomains[1].Solution));
                    using (sw = File.AppendText(fileNameIterations))
                        sw.WriteLine(iterationCount[i - simulationStartFrom].ToString());
                    using (sw = File.CreateText(fileNameTimes))
                    {
                        sw.WriteLine(String.Format("Elements: {0}", times["element"].ToString()));
                        sw.WriteLine(String.Format("Factorize: {0}", times["factorize"].ToString()));
                        sw.WriteLine(String.Format("Solution: {0}", times["solution"].ToString()));
                        sw.WriteLine(String.Format("All: {0}", times["all"].ToString()));
                        sw.WriteLine(String.Format("Total iterations: {0}", totalIterations.ToString()));
                    }

                    foreach (var d in model.NodalDOFsDictionary)
                    {
                        if (d.Value.ContainsKey(DOFType.X) && d.Value[DOFType.X] > -1) dofValuesPerPhi[d.Key][0][j] = subdomains[1].Solution[d.Value[DOFType.X]];
                        if (d.Value.ContainsKey(DOFType.Y) && d.Value[DOFType.Y] > -1) dofValuesPerPhi[d.Key][1][j] = subdomains[1].Solution[d.Value[DOFType.Y]];
                        if (d.Value.ContainsKey(DOFType.Z) && d.Value[DOFType.Z] > -1) dofValuesPerPhi[d.Key][2][j] = subdomains[1].Solution[d.Value[DOFType.Z]];
                    }

                    currentSimulation++;
                }
                foreach (var d in model.NodalDOFsDictionary)
                {
                    if (d.Value.ContainsKey(DOFType.X) && d.Value[DOFType.X] > -1) vrfPerNPtsVRF[d.Key][0][i] = StandardDeviation(dofValuesPerPhi[d.Key][0]) / coefficientsProvider.SpectrumStandardDeviation;
                    if (d.Value.ContainsKey(DOFType.Y) && d.Value[DOFType.Y] > -1) vrfPerNPtsVRF[d.Key][1][i] = StandardDeviation(dofValuesPerPhi[d.Key][1]) / coefficientsProvider.SpectrumStandardDeviation;
                    if (d.Value.ContainsKey(DOFType.Z) && d.Value[DOFType.Z] > -1) vrfPerNPtsVRF[d.Key][2][i] = StandardDeviation(dofValuesPerPhi[d.Key][2]) / coefficientsProvider.SpectrumStandardDeviation;
                }
            }
            times["all"] = DateTime.Now - start;

            vrfs = CalculateVRF(vrfPerNPtsVRF);
            //File.WriteAllLines(String.Format(@"{0}-{1}.txt", fileNameForLogging, expansionOrder), values);
        }

        private Dictionary<int, double[]> CalculateVRF(Dictionary<int, double[][]> vrfPerNPts)
        {
            var vrfs = new Dictionary<int, double[]>();
            foreach (var d in model.NodalDOFsDictionary)
                vrfs.Add(d.Key, new double[3]);

            double dw = coefficientsProvider.Wu / (double)coefficientsProvider.FrequencyIntervals;
            foreach (var d in model.NodalDOFsDictionary)
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < coefficientsProvider.NPtsVRF; k++)
                        vrfs[d.Key][j] += vrfPerNPts[d.Key][j][k] * vrfPerNPts[d.Key][j][k] * coefficientsProvider.SffTarget[k] * dw;
                    vrfs[d.Key][j] = Math.Sqrt(vrfs[d.Key][j]);
                }

            return vrfs;
        }

        private double StandardDeviation(double[] values)
        {
            double avg = values.Average();
            double sumOfSquares = values.Sum(x =>
                {
                    double y = x - avg;
                    return y * y;
                });
            return Math.Sqrt(sumOfSquares / (values.Length - 1));
        }

        private string GetCSVText(IVector iVector)
        {
            var s = string.Empty;
            foreach (var val in model.Subdomains[0].GlobalNodalDOFsDictionary)
                s += (val.Value[DOFType.Y] > -1 ? iVector[val.Value[DOFType.Y]].ToString() : "0") + ";";
            return s;
        }

        private void MakePreconditioner(int simulation)
        {
            int matrixNo = matrixOrder[simulation + blockSize / 2];
            string name = String.IsNullOrWhiteSpace(stiffnessMatrixPath) ? "K" : stiffnessMatrixPath;
            string path = Path.GetDirectoryName(name);
            string nameOnly = Path.GetFileNameWithoutExtension(name);
            string ext = Path.GetExtension(name);

            foreach (var sub in subdomains)
            {
                var m = new SkylineMatrix2D(new int[0]);
                m.ReadFromFile(String.Format("{0}\\{1}Sub{3}Sim{4}{2}", path, nameOnly, ext, sub.Key, matrixNo));
                m.Factorize(1e-8, new List<IVector>(), new List<int>());
                if (factorizedMatrices.ContainsKey(sub.Key))
                    factorizedMatrices[sub.Key] = m;
                else
                    factorizedMatrices.Add(sub.Key, m);
            }
        }

        private void WriteMatricesToFile(int simulation)
        {
            //if (stiffnessMatrixProductionMode == StiffnessMatrixProductionMode.Normal || stiffnessMatrixProductionMode == StiffnessMatrixProductionMode.LoadFromDiskAndCalculate) return;
            string name = String.IsNullOrWhiteSpace(stiffnessMatrixPath) ? "K" : stiffnessMatrixPath;

            string path = Path.GetDirectoryName(name);
            string nameOnly = Path.GetFileNameWithoutExtension(name);
            string ext = Path.GetExtension(name);

            foreach (var sub in subdomains)
                ((IFileWriteable)sub.Value.Matrix).WriteToFile(String.Format(@"{0}\{1}Sub{3}Sim{4}{2}", path, nameOnly, ext, sub.Key, simulation));
        }

        #endregion
    }
}
