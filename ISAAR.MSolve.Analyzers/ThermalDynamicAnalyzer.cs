using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using System.Diagnostics;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    /// <summary>
    /// 
    /// Authors: Yannis Kalogeris
    /// </summary>
    public class ThermalDynamicAnalyzer : IAnalyzer, INonLinearParentAnalyzer
    {
        private readonly double beta, timeStep, totalTime;
        private Dictionary<int, Vector> rhs = new Dictionary<int, Vector>();
        private Dictionary<int, Vector> rhsPrevious = new Dictionary<int, Vector>();
        private Dictionary<int, Vector> Temperature = new Dictionary<int, Vector>();
        private Dictionary<int, Vector> uu = new Dictionary<int, Vector>();
        private Dictionary<int, Vector> CapacityTimesTemperature = new Dictionary<int, Vector>();
        private Dictionary<int, Vector> uc = new Dictionary<int, Vector>();
        private Dictionary<int, Vector> ConductivityTimesTemperature = new Dictionary<int, Vector>();

        //private readonly Dictionary<int, IImplicitIntegrationAnalyzerLog> resultStorages =
        //    new Dictionary<int, IImplicitIntegrationAnalyzerLog>();
        private readonly Dictionary<int, ImplicitIntegrationAnalyzerLog> resultStorages =
            new Dictionary<int, ImplicitIntegrationAnalyzerLog>();
        private readonly IDictionary<int, ILinearSystem> subdomains;
        private readonly IImplicitIntegrationProvider provider;
        private IAnalyzer childAnalyzer;
        private IAnalyzer parentAnalyzer;

        public ThermalDynamicAnalyzer(IImplicitIntegrationProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ILinearSystem> subdomains,
            double beta, double timeStep, double totalTime)
        {
            this.provider = provider;
            this.childAnalyzer = embeddedAnalyzer;
            this.subdomains = subdomains;
            this.beta = beta;
            this.timeStep = timeStep;
            this.totalTime = totalTime;
            this.childAnalyzer.ParentAnalyzer = this;
        }

        public Dictionary<int, ImplicitIntegrationAnalyzerLog> ResultStorages { get { return resultStorages; } }


        public void ResetSolutionVectors()
        {
            foreach (ILinearSystem subdomain in subdomains.Values)
                subdomain.Solution.Clear();
        }

        private void InitializeInternalVectors()
        {
            Temperature.Clear();
            rhs.Clear();
            rhsPrevious.Clear();
            uu.Clear();
            CapacityTimesTemperature.Clear();
            uc.Clear();
            ConductivityTimesTemperature.Clear();

            foreach (ILinearSystem subdomain in subdomains.Values)
            {
                int dofs = subdomain.RHS.Length;
                uu.Add(subdomain.ID, new Vector(dofs));
                CapacityTimesTemperature.Add(subdomain.ID, new Vector(dofs));
                uc.Add(subdomain.ID, new Vector(dofs));
                ConductivityTimesTemperature.Add(subdomain.ID, new Vector(dofs));
                Temperature.Add(subdomain.ID, new Vector(dofs));
                rhs.Add(subdomain.ID, new Vector(dofs));
                rhsPrevious.Add(subdomain.ID, new Vector(dofs));

                // Account for initial conditions coming from a previous solution
                subdomain.Solution.CopyTo(Temperature[subdomain.ID].Data, 0);
            }
        }

        private void InitializeMatrices()
        {
            ImplicitIntegrationCoefficients coeffs = new ImplicitIntegrationCoefficients
            {
                Mass = 1 / timeStep ,
                Stiffness = beta
            };
            foreach (ILinearSystem subdomain in subdomains.Values) provider.CalculateEffectiveMatrix(subdomain, coeffs);
        }

        private void InitializeRHSs()
        {
            ImplicitIntegrationCoefficients coeffs = new ImplicitIntegrationCoefficients
            {
                Mass = 0,  //never used
                Stiffness = 0
            };
            foreach (ILinearSystem subdomain in subdomains.Values)
            {
                provider.ProcessRHS(subdomain, coeffs);
                int dofs = subdomain.RHS.Length;
                //rhsPrevious[subdomain.ID] = rhs[subdomain.ID];
                rhs[subdomain.ID] = new Vector(rhs[subdomain.ID].Length);
                rhs[subdomain.ID].CopyFrom(0, subdomain.RHS.Length, subdomain.RHS, 0);
                //for (int i = 0; i < dofs; i++) rhs[subdomain.ID][i] = subdomain.RHS[i];
            }
        }

        private void UpdateResultStorages(DateTime start, DateTime end)
        {
            if (childAnalyzer == null) return;
            foreach (ILinearSystem subdomain in subdomains.Values)
                if (resultStorages.ContainsKey(subdomain.ID))
                    if (resultStorages[subdomain.ID] != null)
                        foreach (var l in childAnalyzer.Logs[subdomain.ID])
                            resultStorages[subdomain.ID].StoreResults(start, end, l);
        }

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog[]> Logs { get { return null; } }
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

        public void Initialize()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Static analyzer must contain an embedded analyzer.");

            InitializeInternalVectors();
            InitializeRHSs();
            childAnalyzer.Initialize();
        }

        private void CalculateRHSImplicit(ILinearSystem subdomain, double[] rhsResult, bool addRHS)
        {
            int id = subdomain.ID;
            provider.MassMatrixVectorProduct(subdomain, Temperature[id], CapacityTimesTemperature[id].Data);
            provider.DampingMatrixVectorProduct(subdomain, Temperature[id], ConductivityTimesTemperature[id].Data);
            if (addRHS) /// what is the meaning of this?
            {
                for (int i = 0; i < subdomain.RHS.Length; i++)
                {
                    rhsResult[i] = (1 - beta) * rhsPrevious[id].Data[i] + beta * rhs[id].Data[i] + 1 / timeStep * CapacityTimesTemperature[id].Data[i] - (1 - beta)* ConductivityTimesTemperature[id].Data[i];
                }
            }
            else
            {
                for (int i = 0; i < subdomain.RHS.Length; i++)
                {
                    rhsResult[i] = (1 - beta) * rhsPrevious[id].Data[i] + beta * rhs[id].Data[i] + 1 / timeStep * CapacityTimesTemperature[id].Data[i] - (1 - beta)  * ConductivityTimesTemperature[id].Data[i];
                }
            }
            rhsPrevious[subdomain.ID] = rhs[subdomain.ID];
        }

        private void CalculateRHSImplicit()
        {
            foreach (ILinearSystem subdomain in subdomains.Values)
            {
                CalculateRHSImplicit(subdomain, ((Vector)subdomain.RHS).Data, true);
            }
        }

        public void Solve()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Static analyzer must contain an embedded analyzer.");

            int numTimeSteps = (int)(totalTime / timeStep);
            for (int t = 0; t < numTimeSteps; ++t)
            {
                Debug.WriteLine("Newmark step: {0}", t);
                provider.GetRHSFromHistoryLoad(t);
                InitializeRHSs();
                // ProcessRHS
                CalculateRHSImplicit();
                DateTime start = DateTime.Now;
                childAnalyzer.Solve();
                DateTime end = DateTime.Now;
                UpdateTemperature(t);
                UpdateResultStorages(start, end);
            }

            //childAnalyzer.Solve();
        }

        private void UpdateTemperature(int timeStep)
        {
            foreach (ILinearSystem subdomain in subdomains.Values)
            {
                int id = subdomain.ID;
                subdomain.Solution.CopyTo(Temperature[id].Data, 0);

            }
        }

        public void BuildMatrices()
        {
             InitializeMatrices();
        }
        #endregion

        #region INonLinearParentAnalyzer Members

        public double[] GetOtherRHSComponents(ILinearSystem subdomain, IVector currentSolution)
        {
 
            Vector result = new Vector(subdomain.Solution.Length);
            Vector temp = new Vector(subdomain.Solution.Length);
            Vector tempResult = new Vector(subdomain.Solution.Length);

            currentSolution.CopyTo(temp.Data, 0);
            temp.Scale(1);
            provider.MassMatrixVectorProduct(subdomain, temp, tempResult.Data);
            result.Add(tempResult);

            return result.Data;
        }

        #endregion
    }
}

