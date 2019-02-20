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
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;
using Troschuetz.Random.Distributions.Continuous;

namespace ISAAR.MSolve.Analyzers
{
    public class VRFMonteCarloAnalyzerWithStochasticMaterial_v2 : IParentAnalyzer
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
        //private readonly StiffnessMatrixProductionMode stiffnessMatrixProductionMode = StiffnessMatrixProductionMode.Normal;
        private readonly string stiffnessMatrixPath = String.Empty;
        private readonly string randomsReadFileName = String.Empty;
        private readonly string fileNameForLogging = "monteCarlo";
        private readonly PowerSpectrumTargetEvaluatorCoefficientsProvider coefficientsProvider;
        private readonly IStochasticMaterialCoefficientsProvider approximateCoefficientsProvider; 
        private readonly List<int> matrixOrder = new List<int>();
        private readonly List<double> matrixMagnitudes = new List<double>();
        private Dictionary<int, double[]> vrfs;

        public IDictionary<int, LdlSkyline> FactorizedMatrices { get; } = new Dictionary<int, LdlSkyline>();

        public VRFMonteCarloAnalyzerWithStochasticMaterial_v2(Model_v2 model, IAnalyzerProvider_v2 provider, 
            IChildAnalyzer embeddedAnalyzer, ISolver_v2 solver, 
            PowerSpectrumTargetEvaluatorCoefficientsProvider coefficientsProvider, 
            IStochasticMaterialCoefficientsProvider approximateCoefficientsProvider, int expansionOrder, int simulations)
        {
            this.ChildAnalyzer = embeddedAnalyzer;
            this.provider = provider;
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.expansionOrder = expansionOrder;
            this.ChildAnalyzer.ParentAnalyzer = this;
            //this.matrices = new Dictionary<int, IMatrix<double>>(subdomains.Count);
            this.matrices = new Dictionary<int, IMatrix>[expansionOrder + 1];
            this.coefficientsProvider = coefficientsProvider;
            this.approximateCoefficientsProvider = approximateCoefficientsProvider;
            this.simulations = coefficientsProvider.NPhi * coefficientsProvider.NPtsVRF;
            //this.stochasticDomain = stochasticDomain;
        }

        public VRFMonteCarloAnalyzerWithStochasticMaterial_v2(Model_v2 model, IAnalyzerProvider_v2 provider, 
            IChildAnalyzer embeddedAnalyzer, ISolver_v2 solver, 
            PowerSpectrumTargetEvaluatorCoefficientsProvider coefficientsProvider, 
            IStochasticMaterialCoefficientsProvider approximateCoefficientsProvider,
            int expansionOrder, int simulations, string fileNameForLogging)
            : this(model, provider, embeddedAnalyzer, solver, coefficientsProvider, approximateCoefficientsProvider, 
                  expansionOrder, simulations)
        {
            this.fileNameForLogging = fileNameForLogging;
        }

        public VRFMonteCarloAnalyzerWithStochasticMaterial_v2(Model_v2 model, IAnalyzerProvider_v2 provider, 
            IChildAnalyzer embeddedAnalyzer, ISolver_v2 solver, 
            PowerSpectrumTargetEvaluatorCoefficientsProvider coefficientsProvider, 
            IStochasticMaterialCoefficientsProvider approximateCoefficientsProvider,
            int expansionOrder, int simulations, 
            StiffnessMatrixProductionMode stiffnessMatrixProductionMode, string fileNameForLogging, string stiffnessMatrixPath)
            : this(model, provider, embeddedAnalyzer, solver, coefficientsProvider, approximateCoefficientsProvider, 
                  expansionOrder, simulations, fileNameForLogging)
        {
            this.stiffnessMatrixPath = stiffnessMatrixPath;
            //this.stiffnessMatrixProductionMode = stiffnessMatrixProductionMode;
        }

        public VRFMonteCarloAnalyzerWithStochasticMaterial_v2(Model_v2 model, IAnalyzerProvider_v2 provider, 
            IChildAnalyzer embeddedAnalyzer, ISolver_v2 solver, 
            PowerSpectrumTargetEvaluatorCoefficientsProvider coefficientsProvider, 
            IStochasticMaterialCoefficientsProvider approximateCoefficientsProvider,
            int expansionOrder, int simulations, int blockSize, 
            StiffnessMatrixProductionMode stiffnessMatrixProductionMode, string fileNameForLogging, string stiffnessMatrixPath)
            : this(model, provider, embeddedAnalyzer, solver, coefficientsProvider, approximateCoefficientsProvider, 
                  expansionOrder, simulations, stiffnessMatrixProductionMode, fileNameForLogging, stiffnessMatrixPath)
        {
            this.blockSize = blockSize;
        }

        public VRFMonteCarloAnalyzerWithStochasticMaterial_v2(Model_v2 model, IAnalyzerProvider_v2 provider, 
            IChildAnalyzer embeddedAnalyzer, ISolver_v2 solver, 
            PowerSpectrumTargetEvaluatorCoefficientsProvider coefficientsProvider, 
            IStochasticMaterialCoefficientsProvider approximateCoefficientsProvider,
            int expansionOrder, int simulations, int blockSize, 
            StiffnessMatrixProductionMode stiffnessMatrixProductionMode, string fileNameForLogging, string stiffnessMatrixPath, 
            string randomsReadFileName, int simulationStartFrom)
            : this(model, provider, embeddedAnalyzer, solver, coefficientsProvider, approximateCoefficientsProvider, 
                  expansionOrder, simulations, blockSize, stiffnessMatrixProductionMode, fileNameForLogging, stiffnessMatrixPath)
        {
            this.randomsReadFileName = randomsReadFileName;
            this.simulationStartFrom = simulationStartFrom;
        }

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog_v2[]> Logs { get; } = new Dictionary<int, IAnalyzerLog_v2[]>();

        public IChildAnalyzer ChildAnalyzer { get; set; }

        public void BuildMatrices()
        {
            if (ChildAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");
            provider.Reset();
            ChildAnalyzer.BuildMatrices();
        }

        public void Initialize(bool isFirstAnalysis)
        {
            if (ChildAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");

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
            if (ChildAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");

            SolveNormal();
        }

        private void SolveNormal()
        {
            //int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[150][DOFType.Y];
            //int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[84][DOFType.Y];
            int dofNo = model.Subdomains[0].FreeDofOrdering.FreeDofs[model.Nodes[80], DOFType.X];
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
            foreach (var d in model.Nodes.Select(node => node.ID))
                dofValuesPerPhi.Add(d, new double[][] { new double[coefficientsProvider.NPhi], new double[coefficientsProvider.NPhi], new double[coefficientsProvider.NPhi] });
            var vrfPerNPtsVRF = new Dictionary<int, double[][]>(model.Nodes.Count);
            foreach (var d in model.Nodes.Select(node => node.ID))
                vrfPerNPtsVRF.Add(d, new double[][] { new double[coefficientsProvider.NPtsVRF], new double[coefficientsProvider.NPtsVRF], new double[coefficientsProvider.NPtsVRF] });

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
                    ChildAnalyzer.Initialize(false);
                    times["factorize"] += DateTime.Now - e;
                    GCSettings.LatencyMode = GCLatencyMode.LowLatency;
                    e = DateTime.Now;
                    ChildAnalyzer.Solve();

                    times["solution"] += DateTime.Now - e;
                    GCSettings.LatencyMode = GCLatencyMode.Batch;
                    values[i - simulationStartFrom] = linearSystems[1].Solution[dofNo].ToString();
                    //values[i] = subdomains[36].Solution[dofNo].ToString();

                    //values[i] = matrixMagnitudes[i].ToString();
                    using (sw = File.AppendText(fileName))
                        //sw.WriteLine(values[i - simulationStartFrom]);
                        sw.WriteLine(GetCSVText(linearSystems[1].Solution));
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

                    foreach ((INode node, DOFType dofType, int dofIdx) in model.GlobalDofOrdering.GlobalFreeDofs)
                    {
                        double displacement = linearSystems[1].Solution[dofIdx];
                        if (dofType == DOFType.X) dofValuesPerPhi[node.ID][0][j] = displacement;
                        if (dofType == DOFType.Y) dofValuesPerPhi[node.ID][1][j] = displacement;
                        if (dofType == DOFType.Z) dofValuesPerPhi[node.ID][2][j] = displacement;
                    }

                    currentSimulation++;
                }
                foreach ((INode node, DOFType dofType, int globalDofIdx) in model.GlobalDofOrdering.GlobalFreeDofs)
                {
                    int localDofIndex = -1;
                    if (dofType == DOFType.X) localDofIndex = 0;
                    if (dofType == DOFType.Y) localDofIndex = 1;
                    if (dofType == DOFType.Z) localDofIndex = 2;
                    if (localDofIndex != -1)
                    {
                        vrfPerNPtsVRF[node.ID][localDofIndex][i] = StandardDeviation(dofValuesPerPhi[node.ID][localDofIndex]) 
                            / coefficientsProvider.SpectrumStandardDeviation;
                    }
                }
            }
            times["all"] = DateTime.Now - start;

            vrfs = CalculateVRF(vrfPerNPtsVRF);
            //File.WriteAllLines(String.Format(@"{0}-{1}.txt", fileNameForLogging, expansionOrder), values);
        }

        private Dictionary<int, double[]> CalculateVRF(Dictionary<int, double[][]> vrfPerNPts)
        {
            var vrfs = new Dictionary<int, double[]>();
            foreach (var d in model.Nodes.Select(node => node.ID))
                vrfs.Add(d, new double[3]);

            double dw = coefficientsProvider.Wu / (double)coefficientsProvider.FrequencyIntervals;
            foreach (var d in model.Nodes.Select(node => node.ID))
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < coefficientsProvider.NPtsVRF; k++)
                        vrfs[d][j] += vrfPerNPts[d][j][k] * vrfPerNPts[d][j][k] * coefficientsProvider.SffTarget[k] * dw;
                    vrfs[d][j] = Math.Sqrt(vrfs[d][j]);
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

        private string GetCSVText(IVectorView iVector)
        {
            var s = string.Empty;
            ISubdomain_v2 subdomain = model.Subdomains[0];
            foreach (var node in subdomain.Nodes)
            {
                bool isFree = subdomain.FreeDofOrdering.FreeDofs.TryGetValue(node, DOFType.Y, out int dofIdx);
                s += (isFree ? dofIdx.ToString() : "0") + ";";
            }
            return s;
        }

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

        private void WriteMatricesToFile(int simulation)
        {
            //if (stiffnessMatrixProductionMode == StiffnessMatrixProductionMode.Normal || stiffnessMatrixProductionMode == StiffnessMatrixProductionMode.LoadFromDiskAndCalculate) return;
            string name = String.IsNullOrWhiteSpace(stiffnessMatrixPath) ? "K" : stiffnessMatrixPath;

            string path = Path.GetDirectoryName(name);
            string nameOnly = Path.GetFileNameWithoutExtension(name);
            string ext = Path.GetExtension(name);

            var writer = new RawArraysWriter();
            foreach (var linearSystem in linearSystems) 
            {
                writer.WriteToFile((ISparseMatrix)linearSystem.Value.Matrix, 
                    String.Format(@"{0}\{1}Sub{3}Sim{4}{2}", path, nameOnly, ext, linearSystem.Key, simulation));
            }
        }

        #endregion
    }
}
