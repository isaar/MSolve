using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    public class NewtonRaphsonNonLinearAnalyzer_v2 : IAnalyzer_v2
    {
        private readonly IReadOnlyList<ILinearSystem_v2> linearSystems;
        private readonly INonLinearSubdomainUpdater_v2[] subdomainUpdaters;
        private readonly int increments;
        private int maxIterations = 1000;
        private int numIterationsForMatrixRebuild = 1;
        private readonly double tolerance = 1e-8;
        private double rhsNorm;
        private INonLinearParentAnalyzer_v2 parentAnalyzer = null;
        private readonly IStructuralModel_v2 model;
        private readonly ISolver_v2 solver;
        private readonly INonLinearProvider_v2 provider;
        private readonly Dictionary<int, IVector> rhs = new Dictionary<int, IVector>();
        private readonly Dictionary<int, IVector> u = new Dictionary<int, IVector>();
        private readonly Dictionary<int, IVector> du = new Dictionary<int, IVector>();
        private readonly Dictionary<int, IVector> uPlusdu = new Dictionary<int, IVector>();
        private Vector globalRHS; //TODO: This was originally readonly 
        private readonly Dictionary<int, LinearAnalyzerLogFactory> logFactories = new Dictionary<int, LinearAnalyzerLogFactory>();
        private readonly Dictionary<int, IAnalyzerLog[]> logs = new Dictionary<int, IAnalyzerLog[]>();

        public NewtonRaphsonNonLinearAnalyzer_v2(IStructuralModel_v2 model, ISolver_v2 solver,
            INonLinearSubdomainUpdater_v2[] subdomainUpdaters, /*ISubdomainGlobalMapping_v2[] mappings,*/
            INonLinearProvider_v2 provider, int increments)
        {
            this.model = model;
            this.solver = solver;
            this.subdomainUpdaters = subdomainUpdaters;
            this.linearSystems = solver.LinearSystems;
            this.provider = provider;
            this.increments = increments;

            //InitializeInternalVectors(); // I transfered that to Initialize()
        }

        public int MaxIterations
        {
            set
            {
                if (value > 0) maxIterations = value;
                else throw new ArgumentException("Iterations number cannot be negative or zero");
            }
        }


        public int NumIterationsForMatrixRebuild
        {
            set
            {
                if (value > 0) numIterationsForMatrixRebuild = value;
                else throw new Exception("Iterations number for matrix rebuild cannot be negative or zero");
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
                    l.StoreResults(start, end, u[id].ToLegacyVector());
        }

        public IncrementalDisplacementsLog IncrementalDisplacementsLog { get; set; }
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
            globalRHS = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
            rhs.Clear();
            u.Clear();
            du.Clear();
            uPlusdu.Clear();

            foreach (ILinearSystem_v2 linearSystem in linearSystems)
            {
                int id = linearSystem.Subdomain.ID;
                //int numSudomainDofs = linearSystem.RhsVector.Length;
                IVector r = linearSystem.RhsVector.Copy();
                r.ScaleIntoThis(1 / (double)increments);
                rhs.Add(id, r);
                u.Add(id, linearSystem.CreateZeroVector());
                du.Add(id, linearSystem.CreateZeroVector());
                uPlusdu.Add(id, linearSystem.CreateZeroVector());
                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRHS);
            }
            rhsNorm = provider.RHSNorm(globalRHS);
        }

        private void UpdateInternalVectors()//TODOMaria this is where I should add the calculation of the internal nodal force vector
        {
            globalRHS.Clear(); //TODO: Is it necessary to clear it? It will be overwritten by subdomain vectors in this method
            foreach (ILinearSystem_v2 linearSystem in linearSystems)
            {
                int id = linearSystem.Subdomain.ID;

                //TODO: directly copy into linearSystem.RhsVector and then scale that.
                IVector r = linearSystem.RhsVector.Copy();
                r.ScaleIntoThis(1 / (double)increments);
                rhs[id] = r;
                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRHS);
            }
            rhsNorm = provider.RHSNorm(globalRHS);
        }

        public void Initialize()
        {
            InitializeInternalVectors();
            solver.Initialize();
        }

        private void UpdateRHS(int step)
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems)
            {
                linearSystem.RhsVector.CopyFrom(rhs[linearSystem.Subdomain.ID]);
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
                for (step = 0; step < maxIterations; step++)
                {
                    solver.Solve();
                    errorNorm = rhsNorm != 0 ? CalculateInternalRHS(increment, step) / rhsNorm : 0;// (rhsNorm*increment/increments) : 0;//TODOMaria this calculates the internal force vector and subtracts it from the external one (calculates the residual)
                    Console.WriteLine($"Increment {increment}, iteration {step}: norm2(error) = {errorNorm}");

                    if (step == 0) firstError = errorNorm;

                    if (IncrementalDisplacementsLog != null) IncrementalDisplacementsLog.StoreDisplacements_v2(uPlusdu); 

                    if (errorNorm < tolerance) break;

                    SplitResidualForcesToSubdomains();//TODOMaria scatter residuals to subdomains
                    if ((step + 1) % numIterationsForMatrixRebuild == 0)
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
            //CopySolutionToSubdomains();//TODOMaria Copy current displacement to subdomains
            //            ClearMaterialStresses();
            DateTime end = DateTime.Now;

            // TODO: Logging should be done at each iteration. And it should be done using pull observers
            StoreLogResults(start, end);
        }

        private double CalculateInternalRHS(int currentIncrement, int step)
        {
            globalRHS.Clear(); //TODO: Is it necessary to clear it? It will be overwritten by subdomain vectors in this method
            foreach (ILinearSystem_v2 linearSystem in linearSystems)
            {
                int id = linearSystem.Subdomain.ID;

                if (currentIncrement == 0 && step == 0)
                {
                    //TODO: instead of Clear() and then AddIntoThis(), use only CopyFromVector()
                    du[id].Clear();
                    uPlusdu[id].Clear();
                    du[id].AddIntoThis(linearSystem.Solution);
                    uPlusdu[id].AddIntoThis(linearSystem.Solution);
                    du[id].SubtractIntoThis(u[id]);
                }
                else
                {

                    du[id].AddIntoThis(linearSystem.Solution);
                    //TODO: instead of Clear() and then AddIntoThis(), use only CopyFromVector()
                    uPlusdu[id].Clear();
                    uPlusdu[id].AddIntoThis(u[id]);
                    uPlusdu[id].AddIntoThis(du[id]);
                }
                //Vector<double> internalRHS = (Vector<double>)subdomain.GetRHSFromSolution(u[subdomain.ID], du[subdomain.ID]);
                int subdomainIdx = FindSubdomainIdx(linearSystems, linearSystem);

                //TODO: remove cast
                IVector internalRHS = subdomainUpdaters[subdomainIdx].GetRHSFromSolution(uPlusdu[id], du[id]);//TODOMaria this calculates the internal forces
                provider.ProcessInternalRHS(linearSystem, internalRHS, uPlusdu[id]);//TODOMaria this does nothing
                //(new Vector<double>(u[subdomain.ID] + du[subdomain.ID])).Data);

                if (parentAnalyzer != null)
                {
                    IVector otherRhsComponents = parentAnalyzer.GetOtherRHSComponents(linearSystem, uPlusdu[id]);
                    internalRHS.AddIntoThis(otherRhsComponents);//TODOMaria this does nothing for the static problem
                }

                //new Vector<double>(u[subdomain.ID] + du[subdomain.ID]))));

                linearSystem.RhsVector.Clear(); //TODO: we can copy rhs[subdomain.ID] and then scale it instead of clearing and adding.

                //TODO: the next line adds a vector to itself many times. This is called multiplication and is much faster.
                for (int j = 0; j <= currentIncrement; j++) linearSystem.RhsVector.AddIntoThis(rhs[id]);//TODOMaria this adds the external forces 
                linearSystem.RhsVector.SubtractIntoThis(internalRHS);

                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRHS);
            }
            return provider.RHSNorm(globalRHS);
        }

        private void ClearIncrementalSolutionVector()
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems)
                du[linearSystem.Subdomain.ID].Clear();
        }

        private void SplitResidualForcesToSubdomains()
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems)
            {
                int id = linearSystem.Subdomain.ID;
                linearSystem.RhsVector.Clear(); //TODO: why clear it if it is going to be overwritten immediately afterwards?
                model.GlobalDofOrdering.ExtractVectorSubdomainFromGlobal(linearSystem.Subdomain, globalRHS, 
                    linearSystem.RhsVector);
            }
        }

        private void SaveMaterialStateAndUpdateSolution()
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems)
            {
                int id = linearSystem.Subdomain.ID;
                int subdomainIdx = FindSubdomainIdx(linearSystems, linearSystem);
                subdomainUpdaters[subdomainIdx].UpdateState();
                u[id].AddIntoThis(du[id]);
            }
        }

        private int FindSubdomainIdx(IReadOnlyList<ILinearSystem_v2> allLinearSystems, ILinearSystem_v2 wantedlinearSystem)
        {
            return linearSystems.Select((v, i) => 
                new { System = v, Index = i }).First(x => x.System.Subdomain.ID == wantedlinearSystem.Subdomain.ID).Index;
        }

        //TODO: Remove this method. Analyzers should not mess with the solution vector.
        //This method's purpose is to write the final u vector to the linearSystem.Solution, so that StoreLogResults() can 
        //write it to the loggers. It would be faster and more clear to have StoreLogResults() directly access u.
        //private void CopySolutionToSubdomains()
        //{
        //    foreach (ILinearSystem_v2 subdomain in linearSystems)
        //    {
        //        subdomain.Solution.CopyFrom(u[subdomain.ID]);
        //    }
        //}

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
