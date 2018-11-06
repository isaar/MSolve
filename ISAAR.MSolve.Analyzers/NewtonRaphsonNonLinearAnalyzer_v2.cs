using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    public class NewtonRaphsonNonLinearAnalyzer_v2 : IAnalyzer_v2
    {
        private readonly ILinearSystem_v2[] linearSystems;
        private readonly INonLinearSubdomainUpdater_v2[] subdomainUpdaters;
        private readonly ISubdomainGlobalMapping_v2[] mappings;
        private readonly int increments;
        private readonly int totalDOFs;
        private int maxSteps = 1000;
        private int stepsForMatrixRebuild = 1;
        private readonly double tolerance = 1e-8;
        private double rhsNorm;
        private INonLinearParentAnalyzer_v2 parentAnalyzer = null;
        private readonly ISolver_v2 solver;
        private readonly INonLinearProvider_v2 provider;
        private readonly Dictionary<int, IVector> rhs = new Dictionary<int, IVector>();
        private readonly Dictionary<int, IVector> u = new Dictionary<int, IVector>();
        private readonly Dictionary<int, IVector> du = new Dictionary<int, IVector>();
        private readonly Dictionary<int, IVector> uPlusdu = new Dictionary<int, IVector>();
        private IVector globalRHS; //TODO: This was originally readonly 
        private readonly Dictionary<int, LinearAnalyzerLogFactory> logFactories = new Dictionary<int, LinearAnalyzerLogFactory>();
        private readonly Dictionary<int, IAnalyzerLog[]> logs = new Dictionary<int, IAnalyzerLog[]>();

        public NewtonRaphsonNonLinearAnalyzer_v2(ISolver_v2 solver, 
            INonLinearSubdomainUpdater_v2[] subdomainUpdaters, ISubdomainGlobalMapping_v2[] mappings, 
            INonLinearProvider_v2 provider, int increments, int totalDOFs)
        {
            this.solver = solver;
            this.subdomainUpdaters = subdomainUpdaters;
            this.mappings = mappings;
            this.linearSystems = solver.LinearSystems.Values.ToArray();
            this.provider = provider;
            this.increments = increments;
            this.totalDOFs = totalDOFs;

            InitializeInternalVectors();
        }

        public int SetMaxIterations
        {
            set
            {
                if (value > 0) { this.maxSteps = value; }
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
                    l.StoreResults(start, end, linearSystems[id].Solution.ToLegacyVector());
        }

        public Dictionary<int, LinearAnalyzerLogFactory> LogFactories { get { return logFactories; } }

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog[]> Logs { get { return logs; } }

        public IAnalyzer_v2 ParentAnalyzer
        {
            get { return (IAnalyzer_v2)parentAnalyzer; }
            set { parentAnalyzer = (INonLinearParentAnalyzer_v2)value; }
        }

        public IAnalyzer_v2 ChildAnalyzer
        {
            get { return null; }
            set { throw new InvalidOperationException("Newton-Raphson analyzer cannot contain an embedded analyzer."); }
        }

        public void InitializeInternalVectors()//TODOMaria: this is probably where the initial internal nodal vector is calculated
        {
            globalRHS = Vector.CreateZero(totalDOFs); 
            rhs.Clear();
            u.Clear();
            du.Clear();
            uPlusdu.Clear();

            foreach (ILinearSystem_v2 subdomain in linearSystems)
            {
                int numSudomainDofs = subdomain.RhsVector.Length;
                IVector r = subdomain.RhsVector.Copy();
                r.ScaleIntoThis(1 / (double)increments);
                rhs.Add(subdomain.ID, r);
                u.Add(subdomain.ID, Vector.CreateZero(numSudomainDofs));
                du.Add(subdomain.ID, Vector.CreateZero(numSudomainDofs));
                uPlusdu.Add(subdomain.ID, Vector.CreateZero(numSudomainDofs));
                int subdomainIdx = linearSystems.Select((v, i) => new { System = v, Index = i }).
                    First(x => x.System.ID == subdomain.ID).Index;
                mappings[subdomainIdx].SubdomainToGlobalVector(subdomain.RhsVector, globalRHS);
            }
            rhsNorm = provider.RHSNorm(globalRHS);
        }

        private void UpdateInternalVectors()//TODOMaria this is where I should add the calculation of the internal nodal force vector
        {
            globalRHS.Clear(); //TODO: Is it necessary to clear it? It will be overwritten by subdomain vectors in this method
            foreach (ILinearSystem_v2 subdomain in linearSystems)
            {
                //TODO: directly copy into subdomain.RhsVector and then scale that.
                IVector r = subdomain.RhsVector.Copy();
                r.ScaleIntoThis(1 / (double)increments);
                rhs[subdomain.ID] = r; 
                int subdomainIdx = linearSystems.Select((v, i) => new { System = v, Index = i }).
                    First(x => x.System.ID == subdomain.ID).Index;
                mappings[subdomainIdx].SubdomainToGlobalVector(subdomain.RhsVector, globalRHS);
            }
            rhsNorm = provider.RHSNorm(globalRHS);
        }

        public void Initialize()
        {
            solver.Initialize();
        }

        private void UpdateRHS(int step)
        {
            foreach (ILinearSystem_v2 subdomain in linearSystems)
            {
                subdomain.RhsVector.CopyFrom(rhs[subdomain.ID]);
                //subdomainRHS.Multiply(step + 1);
            }
        }

        public void Solve()
        {
            InitializeLogs();

            DateTime start = DateTime.Now;
            UpdateInternalVectors();//TODOMaria this divides the externally applied load by the number of increments and scatters it to all subdomains and stores it in the class subdomain dictionary and total external load vector
            for (int increment = 0; increment < increments; increment++)
            {
                double errorNorm = 0;
                ClearIncrementalSolutionVector();//TODOMaria this sets du to 0
                UpdateRHS(increment);//TODOMaria this copies the residuals stored in the class dictionary to the subdomains

                double firstError = 0;
                int step = 0;
                for (step = 0; step < maxSteps; step++)
                {
                    solver.Solve();
                    errorNorm = rhsNorm != 0 ? CalculateInternalRHS(increment, step) / rhsNorm : 0;// (rhsNorm*increment/increments) : 0;//TODOMaria this calculates the internal force vector and subtracts it from the external one (calculates the residual)
                    if (step == 0) firstError = errorNorm;
                    if (errorNorm < tolerance) break;

                    SplitResidualForcesToSubdomains();//TODOMaria scatter residuals to subdomains
                    if ((step + 1) % stepsForMatrixRebuild == 0)
                    {
                        provider.Reset();
                        BuildMatrices();

                        // The next was used to force the SkylineSolver to factorize the matrix. 
                        // The ISolver.Initialize() could potentially perform actions that must not be repeated or are too expesive.
                        //solver.Initialize(); 
                    }
                }
                Debug.WriteLine("NR {0}, first error: {1}, exit error: {2}", step, firstError, errorNorm);
                SaveMaterialStateAndUpdateSolution();
            }
            CopySolutionToSubdomains();//TODOMaria Copy current displacement to subdomains
            //            ClearMaterialStresses();
            DateTime end = DateTime.Now;

            StoreLogResults(start, end);
        }

        private double CalculateInternalRHS(int currentIncrement, int step)
        {
            globalRHS.Clear(); //TODO: Is it necessary to clear it? It will be overwritten by subdomain vectors in this method
            foreach (ILinearSystem_v2 subdomain in linearSystems)
            {
                if (currentIncrement == 0 && step == 0)
                {
                    //TODO: instead of Clear() and then AddIntoThis(), use only CopyFromVector()
                    du[subdomain.ID].Clear(); 
                    uPlusdu[subdomain.ID].Clear();
                    du[subdomain.ID].AddIntoThis(subdomain.Solution);
                    uPlusdu[subdomain.ID].AddIntoThis(subdomain.Solution);
                    du[subdomain.ID].SubtractIntoThis(u[subdomain.ID]);
                }
                else
                {

                    du[subdomain.ID].AddIntoThis(subdomain.Solution);
                    //TODO: instead of Clear() and then AddIntoThis(), use only CopyFromVector()
                    uPlusdu[subdomain.ID].Clear();
                    uPlusdu[subdomain.ID].AddIntoThis(u[subdomain.ID]);
                    uPlusdu[subdomain.ID].AddIntoThis(du[subdomain.ID]);
                }
                //Vector<double> internalRHS = (Vector<double>)subdomain.GetRHSFromSolution(u[subdomain.ID], du[subdomain.ID]);
                int subdomainIdx = linearSystems.Select((v, i) => new { System = v, Index = i }).
                    First(x => x.System.ID == subdomain.ID).Index;

                //TODO: remove cast
                IVector internalRHS = subdomainUpdaters[subdomainIdx].GetRHSFromSolution(uPlusdu[subdomain.ID], du[subdomain.ID]);//TODOMaria this calculates the internal forces
                provider.ProcessInternalRHS(subdomain, internalRHS, uPlusdu[subdomain.ID]);//TODOMaria this does nothing
                //(new Vector<double>(u[subdomain.ID] + du[subdomain.ID])).Data);

                if (parentAnalyzer != null)
                {
                    IVector otherRhsComponents = parentAnalyzer.GetOtherRHSComponents(subdomain, uPlusdu[subdomain.ID]);
                    internalRHS.AddIntoThis((Vector)otherRhsComponents);//TODOMaria this does nothing for the static problem
                }

                //new Vector<double>(u[subdomain.ID] + du[subdomain.ID]))));

                subdomain.RhsVector.Clear(); //TODO: we can copy rhs[subdomain.ID] and then scale it instead of clearing and adding.

                //TODO: the next line adds a vector to itself many times. This is called multiplication and is much faster.
                for (int j = 0; j <= currentIncrement; j++) subdomain.RhsVector.AddIntoThis(rhs[subdomain.ID]);//TODOMaria this adds the external forces 
                subdomain.RhsVector.SubtractIntoThis(internalRHS);

                subdomainIdx = linearSystems.Select((v, i) => new { System = v, Index = i }).
                    First(x => x.System.ID == subdomain.ID).Index; //TODO: this is already computed. Remove it
                mappings[subdomainIdx].SubdomainToGlobalVector(subdomain.RhsVector, globalRHS);
            }
            return provider.RHSNorm(globalRHS); 
        }

        private void ClearIncrementalSolutionVector()
        {
            foreach (ILinearSystem_v2 subdomain in linearSystems)
                du[subdomain.ID].Clear();
        }

        private void SplitResidualForcesToSubdomains()
        {
            foreach (ILinearSystem_v2 subdomain in linearSystems)
            {
                subdomain.RhsVector.Clear(); //TODO: why clear it if it is going to be overwritten immediately afterwards?
                int subdomainIdx = linearSystems.Select((v, i) => new { System = v, Index = i }).
                    First(x => x.System.ID == subdomain.ID).Index;
                mappings[subdomainIdx].SplitGlobalVectorToSubdomain(globalRHS, subdomain.RhsVector);
            }
        }

        private void SaveMaterialStateAndUpdateSolution()
        {
            foreach (ILinearSystem_v2 subdomain in linearSystems)
            {
                int subdomainIdx = linearSystems.Select((v, i) => new { System = v, Index = i }).
                    First(x => x.System.ID == subdomain.ID).Index;
                subdomainUpdaters[subdomainIdx].UpdateState();
                u[subdomain.ID].AddIntoThis(du[subdomain.ID]);
            }
        }

        //TODO: Remove this method. Analyzers should not mess with the solution vector.
        private void CopySolutionToSubdomains()
        {
            foreach (ILinearSystem_v2 subdomain in linearSystems)
            {
                subdomain.Solution.CopyFrom(u[subdomain.ID]);
            }
        }

        //private void ClearMaterialStresses()
        //{
        //    foreach (ILinearSystem subdomain in linearSystems)
        //        subdomainUpdaters[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].ResetState();
        //}

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
