using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Interfaces;
using System.Collections;
using System.Linq;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.Commons;

namespace ISAAR.MSolve.Analyzers
{
    public class DisplacementControlNonLinearAnalyzer_v2 : IAnalyzer_v2
    {
        private readonly IReadOnlyList<ILinearSystem_v2> linearSystems;
        private readonly INonLinearSubdomainUpdater_v2[] subdomainUpdaters;
        private readonly IEquivalentLoadsAssembler_v2[] equivalentLoadsAssemblers;
        private readonly int increments;
        private int maxiterations = 1000;
        private int stepsForMatrixRebuild = 1;
        private readonly double tolerance = 1e-3;
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

        public DisplacementControlNonLinearAnalyzer_v2(IStructuralModel_v2 model, ISolver_v2 solver,
            INonLinearProvider_v2 provider, INonLinearSubdomainUpdater_v2[] subdomainUpdaters,
            IEquivalentLoadsAssembler_v2[] equivalentLoadsAssemblers,  int increments)
        {
            this.model = model;
            this.solver = solver;
            this.provider = provider;
            this.subdomainUpdaters = subdomainUpdaters;
            this.equivalentLoadsAssemblers = equivalentLoadsAssemblers;
            this.linearSystems = solver.LinearSystems;
            this.increments = increments;

            //InitializeInternalVectors();
        }

        public int MaxIterations
        {
            set
            {
                if (value > 0) { this.maxiterations = value; }
                else { throw new Exception("Max iterations number cannot be negative or zero"); }
            }
        }


        public int NumIterationsForMatrixRebuild
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
                    l.StoreResults(start, end, u[id].ToLegacyVector());
        }

        public Dictionary<int, LinearAnalyzerLogFactory> LogFactories => logFactories;

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog[]> Logs => logs;

        public IAnalyzer_v2 ParentAnalyzer
        {
            get => parentAnalyzer;
            set => parentAnalyzer = (INonLinearParentAnalyzer_v2)value; //TODO: remove this cast. Right it only serves as a check
        }

        public IAnalyzer_v2 ChildAnalyzer
        {
            get => null;
            set => throw new InvalidOperationException("Displacement control analyzer cannot contain an embedded analyzer.");
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
                int idx = FindSubdomainIdx(linearSystems, linearSystem);

                IVector r = linearSystem.RhsVector.Copy();
                r.ScaleIntoThis(1 / (double)increments);
                rhs.Add(id, r);
                u.Add(id, linearSystem.CreateZeroVector());
                du.Add(id, linearSystem.CreateZeroVector());
                uPlusdu.Add(id, linearSystem.CreateZeroVector());
                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRHS);
                subdomainUpdaters[idx].ScaleConstraints(1 / (double)increments);
            }
            rhsNorm = provider.RHSNorm(globalRHS);
        }

        private void UpdateInternalVectors()//TODOMaria this is where I should add the calculation of the internal nodal force vector
        {
            globalRHS.Clear();
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

        private void UpdateRHS(int iteration)
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems)
            {
                linearSystem.RhsVector.CopyFrom(rhs[linearSystem.Subdomain.ID]);
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

                        // The next was used to force the SkylineSolver to factorize the matrix. 
                        // The ISolver.Initialize() could potentially perform actions that must not be repeated or are too expesive.
                        //solver.Initialize(); 
                    }
                }
                Debug.WriteLine("NR {0}, first error: {1}, exit error: {2}", iteration, firstError, errorNorm);
                SaveMaterialStateAndUpdateSolution();
            }
            //CopySolutionToSubdomains(); //TODO: this might have been necessary for the user to get output, but it is bad design
            DateTime end = DateTime.Now;
            StoreLogResults(start, end);
        }

        private double CalculateInternalRHS(int currentIncrement, int iteration)
        {
            globalRHS.Clear();
            foreach (ILinearSystem_v2 linearSystem in linearSystems)
            {
                int id = linearSystem.Subdomain.ID;

                if (currentIncrement == 0 && iteration == 0)
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

        private void AddEquivalentNodalLoadsToRHS(int currentIncrement, int iteration)
        {
            if (iteration != 0)
                return;

            foreach (ILinearSystem_v2 linearSystem in linearSystems)
            {
                int id = linearSystem.Subdomain.ID;
                int idx = FindSubdomainIdx(linearSystems, linearSystem);

                double scalingFactor = 1; //((double)currentIncrement + 2) / (currentIncrement + 1); //2; //
                IVector equivalentNodalLoads = equivalentLoadsAssemblers[idx].GetEquivalentNodalLoads(u[id], scalingFactor);
                linearSystem.RhsVector.SubtractIntoThis(equivalentNodalLoads);

                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRHS);
            }
        }

        private void ScaleSubdomainConstraints(int currentIncrement)
        {
            if (currentIncrement == 0)
                return;

            foreach (ILinearSystem_v2 linearSystem in linearSystems)
            {
                int idx = FindSubdomainIdx(linearSystems, linearSystem);
                double scalingFactor = 1; // ((double)currentIncrement + 2) / (currentIncrement + 1);
                subdomainUpdaters[idx].ScaleConstraints(scalingFactor);
            }
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

        private int FindSubdomainIdx(IReadOnlyList<ILinearSystem_v2> allLinearSystems, ILinearSystem_v2 wantedlinearSystem)
        {
            return linearSystems.Select((v, i) =>
                new { System = v, Index = i }).First(x => x.System.Subdomain.ID == wantedlinearSystem.Subdomain.ID).Index;
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