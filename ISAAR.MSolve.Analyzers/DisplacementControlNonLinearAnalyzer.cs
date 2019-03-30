using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using System.Collections;
using System.Linq;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    public class DisplacementControlNonLinearAnalyzer : IAnalyzer
    {
        private readonly ILinearSystem[] linearSystems;
        private readonly INonLinearSubdomainUpdater[] subdomainUpdaters;
        private readonly ISubdomainGlobalMapping[] mappings;
        private readonly IEquivalentLoadsAssembler[] equivalentLoadsAssemblers;
        private readonly int increments;
        private readonly int totalDOFs;
        private int maxiterations = 1000;
        private int stepsForMatrixRebuild = 1;
        private readonly double tolerance = 1e-3;
        private double rhsNorm;
        private INonLinearParentAnalyzer parentAnalyzer = null;
        private readonly ISolver solver;
        private readonly INonLinearProvider provider;
        private readonly Dictionary<int, IVector> rhs = new Dictionary<int, IVector>();
        private readonly Dictionary<int, Vector> u = new Dictionary<int, Vector>();
        private readonly Dictionary<int, Vector> du = new Dictionary<int, Vector>();
        private readonly Dictionary<int, Vector> uPlusdu = new Dictionary<int, Vector>();
        private readonly Vector globalRHS;
        private readonly Dictionary<int, LinearAnalyzerLogFactory> logFactories = new Dictionary<int, LinearAnalyzerLogFactory>();
        private readonly Dictionary<int, IAnalyzerLog[]> logs = new Dictionary<int, IAnalyzerLog[]>();

        public DisplacementControlNonLinearAnalyzer(ISolver solver, ILinearSystem[] linearSystems, INonLinearSubdomainUpdater[] subdomainUpdaters, ISubdomainGlobalMapping[] mappings,
            IEquivalentLoadsAssembler[] equivalentLoadsAssemblers, INonLinearProvider provider, int increments, int totalDOFs)
        {
            this.solver = solver;
            this.subdomainUpdaters = subdomainUpdaters;
            this.mappings = mappings;
            this.equivalentLoadsAssemblers = equivalentLoadsAssemblers;
            this.linearSystems = linearSystems;
            this.provider = provider;
            this.increments = increments;
            this.totalDOFs = totalDOFs;
            this.globalRHS = new Vector(totalDOFs);

            InitializeInternalVectors();
        }

        public int SetMaxIterations
        {
            set
            {
                if (value > 0) { this.maxiterations = value; }
                else { throw new Exception("Iterations number cannot be negative or zero"); }
            }
        }


        public int SetIterationsForMatrixRebuild
        {
            set
            {
                if (value > 0) { this.stepsForMatrixRebuild = value; }
                else { throw new Exception("Iterations number for matrix rebuild cannot be negative or zero"); }
            }
        }


        private void InitializeLogs()
        {
            logs.Clear();
            foreach (int id in logFactories.Keys) logs.Add(id, logFactories[id].CreateLogs());
        }

        private void StoreLogResults(DateTime start, DateTime end)
        {
            foreach (int id in logs.Keys)
                foreach (var l in logs[id])
                    l.StoreResults(start, end, linearSystems[id].Solution);
        }

        public Dictionary<int, LinearAnalyzerLogFactory> LogFactories { get { return logFactories; } }

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog[]> Logs { get { return logs; } }
        public IAnalyzer ParentAnalyzer
        {
            get { return (IAnalyzer)parentAnalyzer; }
            set { parentAnalyzer = (INonLinearParentAnalyzer)value; }
        }

        public IAnalyzer ChildAnalyzer
        {
            get { return null; }
            set { throw new InvalidOperationException("Newton-Raphson analyzer cannot contain an embedded analyzer."); }
        }

        public void InitializeInternalVectors()//TODOMaria: this is probably where the initial internal nodal vector is calculated
        {
            globalRHS.Clear();
            rhs.Clear();
            u.Clear();
            du.Clear();
            uPlusdu.Clear();

            foreach (ILinearSystem subdomain in linearSystems)
            {
                Vector r = new Vector(subdomain.RHS.Length);
                ((Vector)subdomain.RHS).CopyTo(r.Data, 0);
                r.Multiply(1 / (double)increments);
                rhs.Add(subdomain.ID, r);
                u.Add(subdomain.ID, new Vector(subdomain.RHS.Length));
                du.Add(subdomain.ID, new Vector(subdomain.RHS.Length));
                uPlusdu.Add(subdomain.ID, new Vector(subdomain.RHS.Length));
                mappings[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].SubdomainToGlobalVector(((Vector)subdomain.RHS).Data, globalRHS.Data);
                subdomainUpdaters[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].ScaleConstraints(1 / (double)increments);
            }
            rhsNorm = provider.RHSNorm(globalRHS.Data);
        }

        private void UpdateInternalVectors()//TODOMaria this is where I should add the calculation of the internal nodal force vector
        {
            globalRHS.Clear();
            foreach (ILinearSystem subdomain in linearSystems)
            {
                Vector r = new Vector(subdomain.RHS.Length);
                ((Vector)subdomain.RHS).CopyTo(r.Data, 0);
                r.Multiply(1 / (double)increments);
                rhs[subdomain.ID] = r;
                mappings[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].SubdomainToGlobalVector(((Vector)subdomain.RHS).Data, globalRHS.Data);
            }
            rhsNorm = provider.RHSNorm(globalRHS.Data);
        }

        public void Initialize()
        {
            solver.Initialize();
        }

        private void UpdateRHS(int iteration)
        {
            foreach (ILinearSystem subdomain in linearSystems)
            {
                Vector subdomainRHS = ((Vector)subdomain.RHS);
                rhs[subdomain.ID].CopyTo(subdomainRHS.Data, 0);
                //subdomainRHS.Multiply(iteration + 1);
            }
        }

        public void Solve()
        {
            InitializeLogs();

            DateTime start = DateTime.Now;           
            UpdateInternalVectors();
            for (int increment = 0; increment < increments; increment++)
            {
                double errorNorm = 0;
                ClearIncrementalSolutionVector();
                UpdateRHS(increment);
                ScaleSubdomainConstraints(increment);

                double firstError = 0;
                int iteration = 0;
                for (iteration = 0; iteration < maxiterations; iteration++)
                {
                    AddEquivalentNodalLoadsToRHS(increment, iteration);
                    solver.Solve();
                    errorNorm = CalculateInternalRHS(increment, iteration);
                    if (iteration == 0) firstError = errorNorm;
                    if (errorNorm < tolerance) break;

                    SplitResidualForcesToSubdomains();
                    if ((iteration + 1) % stepsForMatrixRebuild == 0)
                    {
                        provider.Reset();
                        BuildMatrices();
                        solver.Initialize();
                    }
                }
                Debug.WriteLine("NR {0}, first error: {1}, exit error: {2}", iteration, firstError, errorNorm);
                SaveMaterialStateAndUpdateSolution();
            }
            CopySolutionToSubdomains();
            DateTime end = DateTime.Now;
            StoreLogResults(start, end);
        }

        private double CalculateInternalRHS(int currentIncrement, int iteration)
        {
            globalRHS.Clear();
            foreach (ILinearSystem subdomain in linearSystems)
            {
                if (currentIncrement == 0 && iteration == 0)
                {
                    Array.Clear(du[subdomain.ID].Data, 0, du[subdomain.ID].Length);
                    Array.Clear(uPlusdu[subdomain.ID].Data, 0, uPlusdu[subdomain.ID].Length);
                    du[subdomain.ID].Add(((Vector)subdomain.Solution));
                    uPlusdu[subdomain.ID].Add(((Vector)subdomain.Solution));
                    du[subdomain.ID].Subtract(u[subdomain.ID]);
                }
                else
                {
                    du[subdomain.ID].Add(((Vector)subdomain.Solution));
                    Array.Clear(uPlusdu[subdomain.ID].Data, 0, uPlusdu[subdomain.ID].Length);
                    uPlusdu[subdomain.ID].Add(u[subdomain.ID]);
                    uPlusdu[subdomain.ID].Add(du[subdomain.ID]);
                }
                //Vector<double> internalRHS = (Vector<double>)subdomain.GetRHSFromSolution(u[subdomain.ID], du[subdomain.ID]);
                Vector internalRHS = (Vector)subdomainUpdaters[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].GetRHSFromSolution(uPlusdu[subdomain.ID], du[subdomain.ID]);//TODOMaria this calculates the internal forces
                provider.ProcessInternalRHS(subdomain, internalRHS.Data, uPlusdu[subdomain.ID].Data);//TODOMaria this does nothing
                //(new Vector<double>(u[subdomain.ID] + du[subdomain.ID])).Data);

                if (parentAnalyzer != null)
                    internalRHS.Add(new Vector(parentAnalyzer.GetOtherRHSComponents(subdomain,
                        uPlusdu[subdomain.ID])));//TODOMaria this does nothing for the static problem
                //new Vector<double>(u[subdomain.ID] + du[subdomain.ID]))));

                Vector subdomainRHS = ((Vector)subdomain.RHS);
                subdomainRHS.Clear();
                for (int j = 0; j <= currentIncrement; j++) subdomainRHS.Add((Vector)rhs[subdomain.ID]);//TODOMaria this adds the external forces 
                subdomainRHS.Subtract(internalRHS);
                mappings[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].SubdomainToGlobalVector(subdomainRHS.Data, globalRHS.Data);
            }
            double providerRHSNorm = provider.RHSNorm(globalRHS.Data);
            return providerRHSNorm;
        }

        private void AddEquivalentNodalLoadsToRHS(int currentIncrement, int iteration)
        {
            if (iteration != 0)
                return;

            foreach (ILinearSystem subdomain in linearSystems)
            {
                Vector subdomainRHS = ((Vector)subdomain.RHS);

                //TODO: Remove the next line. It does nothing.
                var subdomainUpdater = subdomainUpdaters[linearSystems.Select((v, i) => new { System = v, Index = i })
                                                                      .First(x => x.System.ID == subdomain.ID).Index];

                var equivalentLoadsAssembler = equivalentLoadsAssemblers[linearSystems.Select((v, i) => new { System = v, Index = i })
                                                                                      .First(x => x.System.ID == subdomain.ID).Index];

                double scalingFactor = 1; //((double)currentIncrement + 2) / (currentIncrement + 1); //2; //
                var equivalentNodalLoads = (Vector)equivalentLoadsAssembler.GetEquivalentNodalLoads(u[subdomain.ID], scalingFactor);
                subdomainRHS.Subtract(equivalentNodalLoads);

                mappings[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].SubdomainToGlobalVector(subdomainRHS.Data, globalRHS.Data);
            }
        }

        private void ScaleSubdomainConstraints(int currentIncrement)
        {
            if (currentIncrement == 0)
                return;

            foreach (ILinearSystem subdomain in linearSystems)
            {
                Vector subdomainRHS = ((Vector)subdomain.RHS);
                var subdomainUpdater = subdomainUpdaters[linearSystems.Select((v, i) => new { System = v, Index = i })
                                                                      .First(x => x.System.ID == subdomain.ID).Index];

                double scalingFactor = 1; // ((double)currentIncrement + 2) / (currentIncrement + 1);
                subdomainUpdater.ScaleConstraints(scalingFactor);
            }
        }

        private void ClearIncrementalSolutionVector()
        {
            foreach (ILinearSystem subdomain in linearSystems)
                du[subdomain.ID].Clear();
        }

        private void SplitResidualForcesToSubdomains()
        {
            foreach (ILinearSystem subdomain in linearSystems)
            {
                Vector subdomainRHS = ((Vector)subdomain.RHS);
                subdomainRHS.Clear();
                mappings[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].SplitGlobalVectorToSubdomain(globalRHS.Data, subdomainRHS.Data);
            }
        }

        private void SaveMaterialStateAndUpdateSolution()
        {
            foreach (ILinearSystem subdomain in linearSystems)
            {
                subdomainUpdaters[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].UpdateState();
                u[subdomain.ID].Add(du[subdomain.ID]);
            }
        }

        private void CopySolutionToSubdomains()
        {
            foreach (ILinearSystem subdomain in linearSystems)
                u[subdomain.ID].CopyTo(((Vector)subdomain.Solution).Data, 0);
        }

        private void ClearMaterialStresses()
        {
            foreach (ILinearSystem subdomain in linearSystems)
                subdomainUpdaters[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].ResetState();
        }

        public void BuildMatrices()
        {
            if (parentAnalyzer == null)
                throw new InvalidOperationException("This Newton-Raphson non-linear analyzer has no parent.");

            parentAnalyzer.BuildMatrices();
            //solver.Initialize();
        }

        #endregion
    }
}