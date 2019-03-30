using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using System.Collections;
using System.Linq;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Solvers;

namespace ISAAR.MSolve.Analyzers
{
    /// <summary>
    /// Newton Raphson Non-Linear Analyzer for the solution of boundary value problems that arise 
    /// from a FE2 multiscale Analysis and are related to the microstructure model.
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class MicrostructureBvpNRNLAnalyzer : IChildAnalyzer
    {
        private readonly IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems;
        protected readonly IStructuralModel_v2 model;
        private readonly Dictionary<int,NonLinearSubdomainUpdaterWithInitialConditions_v2> subdomainUpdaters;
        //private readonly ISubdomainGlobalMapping[] mappings; 
        private readonly int increments;
        //private readonly int totalDOFs;
        private int maxSteps = 1000;
        private int stepsForMatrixRebuild = 0;
        private readonly double tolerance = 1e-3;
        private double rhsNorm;
        private INonLinearParentAnalyzer_v2 parentAnalyzer = null;
        private readonly ISolver_v2 solver;
        private readonly INonLinearProvider_v2 provider;
        private readonly Dictionary<int, IVector> rhs = new Dictionary<int, IVector>();//comment MS2:apothikevetai se afto h timh (externalLoads/increments) gia kathe subdomain kai apo ekei pernietai opou xreiasthei (p.x. subdomain.RHS)
        private readonly Dictionary<int, IVector> u = new Dictionary<int, IVector>();
        private readonly Dictionary<int, IVector> du = new Dictionary<int, IVector>();
        private readonly Dictionary<int, IVector> uPlusdu = new Dictionary<int, IVector>();
        Dictionary<int, Node_v2> boundaryNodes;
        Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements;
        Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements;
        private readonly Dictionary<int, EquivalentContributionsAssebler_v2> equivalentContributionsAssemblers;
        private Vector globalRhs;
        private readonly Dictionary<int, LinearAnalyzerLogFactory_v2> logFactories = new Dictionary<int, LinearAnalyzerLogFactory_v2>();
        private readonly Dictionary<int, IAnalyzerLog_v2[]> logs = new Dictionary<int, IAnalyzerLog_v2[]>();

        public MicrostructureBvpNRNLAnalyzer(IStructuralModel_v2 model, ISolver_v2 solver, Dictionary<int, NonLinearSubdomainUpdaterWithInitialConditions_v2> subdomainUpdaters,
            INonLinearProvider_v2 provider, int increments, Dictionary<int, IVector> uInitialFreeDOFDisplacementsPerSubdomain,
            Dictionary<int, Node_v2> boundaryNodes, Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements,
            Dictionary<int, EquivalentContributionsAssebler_v2> equivalentContributionsAssemblers)//, ISubdomainGlobalMapping[] mappings)
        {
            this.model = model;
            this.solver = solver;
            this.subdomainUpdaters = subdomainUpdaters;
            //this.mappings = mappings;
            this.linearSystems = solver.LinearSystems;
            this.provider = provider;
            this.increments = increments;
            //this.globalRhs = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs); 
            this.u = uInitialFreeDOFDisplacementsPerSubdomain; //prosthiki MS, TODO:possibly check for compatibility elements format: u.Add(subdomain.ID, new Vector(subdomain.RHS.Length));
                                                               //this.uPlusdu = uInitialFreeDOFDisplacementsPerSubdomain; //prosthiki MS: commented out possible pass by reference
            this.boundaryNodes = boundaryNodes;
            this.initialConvergedBoundaryDisplacements = initialConvergedBoundaryDisplacements;
            this.totalBoundaryDisplacements = totalBoundaryDisplacements;
            this.equivalentContributionsAssemblers = equivalentContributionsAssemblers;

            //InitializeInternalVectors();
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

        //private void StoreLogResults(DateTime start, DateTime end)
        //{
        //    foreach (int id in logs.Keys)
        //        foreach (var l in logs[id])
        //            l.StoreResults(start, end, linearSystems[id].Solution);
        //}

        public Dictionary<int, LinearAnalyzerLogFactory_v2> LogFactories { get { return logFactories; } }

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog_v2[]> Logs { get { return logs; } }

        public IParentAnalyzer ParentAnalyzer// exei diorthothei apo bas v2
        {
            get => parentAnalyzer;
            set => parentAnalyzer = (INonLinearParentAnalyzer_v2)value; //TODO: remove this cast. Now it only serves as a check
        }

        public IAnalyzer ChildAnalyzer
        {
            get { return null; }
            set { throw new InvalidOperationException("Newton-Raphson analyzer cannot contain an embedded analyzer."); }
        }

        public void InitializeInternalVectors()//TODOMaria: this is probably where the initial internal nodal vector is calculated
        {
            globalRhs = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
            rhs.Clear();
            //u.Clear(); prosthiki MS
            du.Clear();
            uPlusdu.Clear(); //prosthiki MS

            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;

                IVector r = linearSystem.RhsVector.Copy();
                r.ScaleIntoThis(1 / (double)increments);
                rhs.Add(id, r); //comment MS: xtizei to rhs, field tou NRNLAnalyzer ok, swsta giati rhs sth dhmiourgia twn linear systems exoume perasei to model.Subdomains[0].Forces
                                //u.Add(subdomain.ID, new Vector(subdomain.RHS.Length)); //prosthiki MS
                du.Add(id, linearSystem.CreateZeroVector());
                var tempcopy = u[id].CopyToArray();
                uPlusdu.Add(id, Vector.CreateFromArray(tempcopy));
                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRhs);

                //Vector tempCopy = new Vector(u[linearSystem.ID].Length);
                //u[linearSystem.ID].Data.CopyTo(tempCopy.Data, 0);
                //uPlusdu.Add(linearSystem.ID, tempCopy); //prosthiki MS
                //mappings[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == linearSystem.ID).Index].SubdomainToGlobalVector(((Vector)linearSystem.RHS).Data, globalRhs.Data);
            }
            rhsNorm = provider.CalculateRhsNorm(globalRhs);
        }

        private void UpdateInternalVectors()//TODOMaria this is where I should add the calculation of the internal nodal force vector
        {
            globalRhs.Clear();
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;

                IVector r = linearSystem.RhsVector.Copy();
                r.ScaleIntoThis(1 / (double)increments);
                rhs[id] = r; //comment MS: xanahtizei to rhs kai to xanaperna sto globalRHS(molis to kane clear opote ok)
                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRhs);
            }
            rhsNorm = provider.CalculateRhsNorm(globalRhs);
        }

        public void Initialize(bool isFirstAnalysis = true)
        {
            InitializeInternalVectors();
            solver.Initialize();
        }

        private void UpdateRHS(int step)
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                linearSystem.RhsVector.CopyFrom(rhs[linearSystem.Subdomain.ID]);                                
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
                UpdateRHS(increment);//comment MS2: apo to rhs[subdomain.ID] pernaei sto subdomain.RHS h fixed timh (externalLoads/increments) (ginetai copy kai oxi add)  AFTO thewreitai RHS sthn prwth iteration
                UpdateRHSForLinearizationContributions(increment + 1, increments);// opote edw pouu uparxei to neo rhs tou kuklou epanalhpsewn prepei na afairountai oi draseis
                double firstError = 0;
                int step = 0;
                for (step = 0; step < maxSteps; step++)
                {
                    solver.Solve();
                    errorNorm = rhsNorm != 0 ? CalculateInternalRHS(increment, step, increments) / rhsNorm : 0;//comment MS2: to subdomain.RHS lamvanei thn timh nIncrement*(externalLoads/increments)-interanalRHS me xrhsh ths fixed timhs apo to rhs[subdomain.ID]
                    if (step == 0) firstError = errorNorm;
                    if (errorNorm < tolerance) break;

                    SplitResidualForcesToSubdomains();//TODOMaria scatter residuals to subdomains
                    if ((step + 1) % stepsForMatrixRebuild == 0)
                    {
                        provider.Reset();
                        BuildMatrices();
                        solver.Initialize();
                    }
                }
                Debug.WriteLine("NR {0}, first error: {1}, exit error: {2}", step, firstError, errorNorm);
                UpdateSolution();
            }
            CopySolutionToSubdomains();//TODOMaria Copy current displacement to subdomains
                                       //            ClearMaterialStresses();
            DateTime end = DateTime.Now;

            //StoreLogResults(start, end);
        }

        private double CalculateInternalRHS(int currentIncrement, int step, int totalIncrements)
        {
            globalRhs.Clear();
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;

                if (currentIncrement == 0 && step == 0)
                {                    
                    du[id].AddIntoThis((linearSystem.Solution));
                    uPlusdu[id].Clear();
                    uPlusdu[id].AddIntoThis(u[id]);
                    uPlusdu[id].AddIntoThis(du[id]);
                }
                else
                {
                    du[id].AddIntoThis((linearSystem.Solution));
                    uPlusdu[id].Clear();
                    uPlusdu[id].AddIntoThis(u[id]);
                    uPlusdu[id].AddIntoThis(du[id]);
                }
                IVector internalRhs = subdomainUpdaters[id].GetRHSFromSolutionWithInitialDisplacemntsEffect(uPlusdu[id], du[id], boundaryNodes,
                initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, currentIncrement + 1, totalIncrements);//TODOMaria this calculates the internal forces
                provider.ProcessInternalRhs(linearSystem.Subdomain, uPlusdu[id], internalRhs);//TODOMaria this does nothing
                                                                                    //(new Vector<double>(u[subdomain.ID] + du[subdomain.ID])).Data);

                if (parentAnalyzer != null)
                {
                    IVector otherRhsComponents = parentAnalyzer.GetOtherRhsComponents(linearSystem, uPlusdu[id]);
                    internalRhs.AddIntoThis(otherRhsComponents);//TODOMaria this does nothing for the static problem
                }//TODOMaria this does nothing for the static problem
                 //new Vector<double>(u[subdomain.ID] + du[subdomain.ID]))));

                linearSystem.RhsVector.Clear();
                for (int j = 0; j <= currentIncrement; j++) linearSystem.RhsVector.AddIntoThis(rhs[id]);//TODOMaria this adds the external forces 
                linearSystem.RhsVector.SubtractIntoThis(internalRhs);
                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRhs);
            }
            double providerRHSNorm = provider.CalculateRhsNorm(globalRhs);
            return providerRHSNorm;
        }

        private void ClearIncrementalSolutionVector()
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                du[id].Clear();
            }
        }

        private void SplitResidualForcesToSubdomains()
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                linearSystem.RhsVector.Clear();
                model.GlobalDofOrdering.ExtractVectorSubdomainFromGlobal(linearSystem.Subdomain, globalRhs, linearSystem.RhsVector);                
            }
        }

        private void UpdateSolution()
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                u[id].AddIntoThis(du[id]);
            }
        }

        public void SaveMaterialState()
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                subdomainUpdaters[id].UpdateState();
            }
        }

        public Dictionary<int, IVector> GetConvergedSolutionVectorsOfFreeDofs()
        {
            return u;
            // return uplusdu einai to idio afou exei ginei molis UpdateSolution
        }

        public Dictionary<int, IVector> GetConvergedIncrementalSolutionVectorsOfFreeDofs()
        {
            return du;
        }


        private void CopySolutionToSubdomains()
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                //TODO:
                //int id = linearSystem.Subdomain.ID;
                //u[id].CopyTo(((Vector)linearSystem.Solution).Data, 0);
            }
        }

        private void ClearMaterialStresses()
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                subdomainUpdaters[id].ResetState();
            }
        }

        public void BuildMatrices()
        {
            if (parentAnalyzer == null)
                throw new InvalidOperationException("This Newton-Raphson non-linear analyzer has no parent.");

            parentAnalyzer.BuildMatrices();
            //solver.Initialize();
        }

        private void UpdateRHSForLinearizationContributions(int nIncrement, int increments)
        {
            globalRhs.Clear();
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                //var equivalentContributionsAssembler = equivalentContributionsAssemblers[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index];
                var equivalentContributionsAssembler = equivalentContributionsAssemblers[id];
                Vector Contribution = equivalentContributionsAssembler.CalculateKfreeprescribedUpMultiplicationForSubdRHSContribution(
                    boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, increments);
                linearSystem.RhsVector.SubtractIntoThis(Contribution);
                //Vector subdomainRHS = ((Vector)linearSystem.RHS);
                //subdomainRHS.Subtract(Contribution);
                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRhs);
                //mappings[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == linearSystem.ID).Index].SubdomainToGlobalVector(subdomainRHS.Data, globalRhs.Data);
            }
            rhsNorm = provider.CalculateRhsNorm(globalRhs);
            //TODOMS: possibly it is not nessesary to update globalRHS (and of course not clear it before the loop) as it is not updated in 177-UpdateRHS(increment) either.
            // and it is used when it is recalculated in CalculateInternalRHS......
        }

        #endregion
    }
}
