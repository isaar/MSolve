using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: It would be clearer and less error prone if incrementing the loading conditions doesn't alter the model's/subdomain's 
//      loads and constraints.
namespace ISAAR.MSolve.Analyzers.NonLinear
{
    public abstract class NonLinearAnalyzerBase : IChildAnalyzer
    {
        //TODO: this should be passed in the constructor by the implementing class and used in Solve().
        //protected readonly double residualTolerance; 

        protected readonly IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems;
        protected readonly int maxIterationsPerIncrement;
        protected readonly IStructuralModel_v2 model;
        protected readonly int numIncrements;
        protected readonly int numIterationsForMatrixRebuild;
        protected readonly INonLinearProvider_v2 provider;
        protected readonly double residualTolerance;
        protected readonly ISolver_v2 solver;
        protected readonly IReadOnlyDictionary<int, INonLinearSubdomainUpdater_v2> subdomainUpdaters;
        protected readonly Dictionary<int, IVector> rhs = new Dictionary<int, IVector>();
        protected readonly Dictionary<int, IVector> u = new Dictionary<int, IVector>();
        protected readonly Dictionary<int, IVector> du = new Dictionary<int, IVector>();
        protected readonly Dictionary<int, IVector> uPlusdu = new Dictionary<int, IVector>();
        protected Vector globalRhs; //TODO: This was originally readonly 
        protected double globalRhsNormInitial; //TODO: This can probably be a local variable.
        protected INonLinearParentAnalyzer_v2 parentAnalyzer = null;

        internal NonLinearAnalyzerBase(IStructuralModel_v2 model, ISolver_v2 solver, INonLinearProvider_v2 provider,
            IReadOnlyDictionary<int, INonLinearSubdomainUpdater_v2> subdomainUpdaters,
            int numIncrements, int maxIterationsPerIncrement, int numIterationsForMatrixRebuild, double residualTolerance)
        {
            this.model = model;
            this.solver = solver;
            this.provider = provider;
            this.subdomainUpdaters = subdomainUpdaters;
            this.linearSystems = solver.LinearSystems;
            this.numIncrements = numIncrements;
            this.maxIterationsPerIncrement = maxIterationsPerIncrement;
            this.numIterationsForMatrixRebuild = numIterationsForMatrixRebuild;
            this.residualTolerance = residualTolerance;
        }

        public Dictionary<int, LinearAnalyzerLogFactory_v2> LogFactories { get; } = new Dictionary<int, LinearAnalyzerLogFactory_v2>();
        public Dictionary<int, IAnalyzerLog_v2[]> Logs { get; } = new Dictionary<int, IAnalyzerLog_v2[]>();

        public TotalDisplacementsPerIterationLog_v2 TotalDisplacementsPerIterationLog { get; set; }
        public Dictionary<int, TotalLoadsDisplacementsPerIncrementLog> IncrementalLogs { get; }
            = new Dictionary<int, TotalLoadsDisplacementsPerIncrementLog>();

        public IParentAnalyzer ParentAnalyzer
        {
            get => parentAnalyzer;
            set => parentAnalyzer = (INonLinearParentAnalyzer_v2)value; //TODO: remove this cast. Now it only serves as a check
        }

        public void BuildMatrices()
        {
            if (parentAnalyzer == null) throw new InvalidOperationException(
                "This Newton-Raphson nonlinear analyzer has no parent.");

            parentAnalyzer.BuildMatrices();
        }

        public void Initialize(bool isFirstAnalysis)
        {
            InitializeInternalVectors();
            solver.Initialize();
        }

        //TODO: Internal Rhs vectors are created and destroyed at each iteration. It would be more efficient to store them as
        //      vectors and then overwrite them.
        protected Dictionary<int, IVector> CalculateInternalRhs(int currentIncrement, int iteration)
        {
            var internalRhsVectors = new Dictionary<int, IVector>();
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
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
                //Vector<double> internalRhs = (Vector<double>)subdomain.GetRhsFromSolution(u[subdomain.ID], du[subdomain.ID]);

                //TODO: remove cast
                IVector internalRhs = subdomainUpdaters[id].GetRhsFromSolution(uPlusdu[id], du[id]);//TODOMaria this calculates the internal forces
                provider.ProcessInternalRhs(linearSystem.Subdomain, uPlusdu[id], internalRhs);//TODOMaria this does nothing
                //(new Vector<double>(u[subdomain.ID] + du[subdomain.ID])).Data);

                if (parentAnalyzer != null)
                {
                    IVector otherRhsComponents = parentAnalyzer.GetOtherRhsComponents(linearSystem, uPlusdu[id]);
                    internalRhs.AddIntoThis(otherRhsComponents);//TODOMaria this does nothing for the static problem
                }

                internalRhsVectors.Add(id, internalRhs);
            }

            return internalRhsVectors; 
        }

        protected double UpdateResidualForcesAndNorm(int currentIncrement, Dictionary<int, IVector> internalRhs)
        {
            globalRhs.Clear();
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;

                linearSystem.RhsVector.Clear(); //TODO: we can copy rhs[subdomain.ID] and then scale it instead of clearing and adding.

                // External forces = loadFactor * total external forces
                //TODO: the next line adds a vector to itself many times. This is called multiplication and is much faster.
                for (int j = 0; j <= currentIncrement; j++) linearSystem.RhsVector.AddIntoThis(rhs[id]);//TODOMaria this adds the external forces 

                // Residual forces = external - internal
                linearSystem.RhsVector.SubtractIntoThis(internalRhs[id]);

                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRhs);
            }
            return provider.CalculateRhsNorm(globalRhs);
        }

        protected void ClearIncrementalSolutionVector()
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values) du[linearSystem.Subdomain.ID].Clear();
        }

        protected virtual void InitializeInternalVectors()//TODOMaria: this is probably where the initial internal nodal vector is calculated
        {
            globalRhs = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
            rhs.Clear();
            u.Clear();
            du.Clear();
            uPlusdu.Clear();

            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;

                IVector r = linearSystem.RhsVector.Copy();
                r.ScaleIntoThis(1 / (double)numIncrements);
                rhs.Add(id, r);
                u.Add(id, linearSystem.CreateZeroVector());
                du.Add(id, linearSystem.CreateZeroVector());
                uPlusdu.Add(id, linearSystem.CreateZeroVector());
                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRhs);
            }
            globalRhsNormInitial = provider.CalculateRhsNorm(globalRhs);
        }

        protected void InitializeLogs()
        {
            Logs.Clear();
            foreach (int id in LogFactories.Keys) Logs.Add(id, LogFactories[id].CreateLogs());
            foreach (var log in IncrementalLogs.Values) log.Initialize();
        }

        protected void SaveMaterialStateAndUpdateSolution()
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                subdomainUpdaters[id].UpdateState();
                u[id].AddIntoThis(du[id]);
            }
        }

        protected void SplitResidualForcesToSubdomains()
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                linearSystem.RhsVector.Clear(); //TODO: why clear it if it is going to be overwritten immediately afterwards?
                model.GlobalDofOrdering.ExtractVectorSubdomainFromGlobal(linearSystem.Subdomain, globalRhs,
                    linearSystem.RhsVector);
            }
        }

        protected void StoreLogResults(DateTime start, DateTime end)
        {
            foreach (int id in Logs.Keys)
                foreach (var l in Logs[id])
                    l.StoreResults(start, end, u[id]);
        }

        protected void UpdateInternalVectors()//TODOMaria this is where I should add the calculation of the internal nodal force vector
        {
            globalRhs.Clear();
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;

                //TODO: directly copy into linearSystem.RhsVector and then scale that.
                IVector r = linearSystem.RhsVector.Copy();
                r.ScaleIntoThis(1 / (double)numIncrements);
                rhs[id] = r;
                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRhs);
            }
            globalRhsNormInitial = provider.CalculateRhsNorm(globalRhs);
        }

        protected void UpdateRhs(int step)
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                linearSystem.RhsVector.CopyFrom(rhs[linearSystem.Subdomain.ID]);
                //linearSystem.RhsVector.Multiply(step + 1);
            }
        }

        //TODO: this should be implemented as a (virtual?) template method. That requires restructuring the helper methods, such 
        //      that only the parts of the algorithm that are different for LoadControl, DisplacementControl, etc, are delegated
        //      to the concrete implementation. E.g. instead of the load incrementing be done by many different helper methods
        //      (e.g. ScaleSubdomainConstraints()), which are called in various parts of the algorithm, we could have an abstract
        //      IncrementLoading() method, where LoadControl would increment the nodal loads, DisplacementControl the 
        //      prescribed displacements, etc. 
        public abstract void Solve();
    }
}
