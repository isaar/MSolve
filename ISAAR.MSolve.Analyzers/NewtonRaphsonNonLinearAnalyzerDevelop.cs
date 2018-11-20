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
using System.Linq; using ISAAR.MSolve.Analyzers.SupportiveClasses;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    public class NewtonRaphsonNonLinearAnalyzerDevelop : IAnalyzer
    {
        private readonly ILinearSystem[] linearSystems;
        private readonly NonLinearSubdomainUpdaterWithInitialConditions[] subdomainUpdaters;
        private readonly ISubdomainGlobalMapping[] mappings;
        private readonly int increments;
        private readonly int totalDOFs;
        private int maxSteps = 1000;
        private int stepsForMatrixRebuild = 0;
        private readonly double tolerance = 1e-3;
        private double rhsNorm;
        private INonLinearParentAnalyzer parentAnalyzer = null;
        private readonly ISolver solver;
        private readonly INonLinearProvider provider;
        private readonly Dictionary<int, IVector> rhs = new Dictionary<int, IVector>();//comment MS2:apothikevetai se afto h timh (externalLoads/increments) gia kathe subdomain kai apo ekei pernietai opou xreiasthei (p.x. subdomain.RHS)
        private readonly Dictionary<int, Vector> u = new Dictionary<int, Vector>(); 
        private readonly Dictionary<int, Vector> du = new Dictionary<int, Vector>();
        private readonly Dictionary<int, Vector> uPlusdu = new Dictionary<int, Vector>();
        Dictionary<int, Node> boundaryNodes;
        Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements;
        Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements;
        private readonly Dictionary<int,EquivalentContributionsAssebler> equivalentContributionsAssemblers;
        private readonly Vector globalRHS;
        private readonly Dictionary<int, LinearAnalyzerLogFactory> logFactories = new Dictionary<int, LinearAnalyzerLogFactory>();
        private readonly Dictionary<int, IAnalyzerLog[]> logs = new Dictionary<int, IAnalyzerLog[]>();

        public NewtonRaphsonNonLinearAnalyzerDevelop(ISolver solver, ILinearSystem[] linearSystems, NonLinearSubdomainUpdaterWithInitialConditions[] subdomainUpdaters, ISubdomainGlobalMapping[] mappings,
            INonLinearProvider provider, int increments, int totalDOFs, Dictionary<int, Vector> uInitialFreeDOFDisplacementsPerSubdomain,
            Dictionary<int, Node> boundaryNodes, Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements,
            Dictionary<int, EquivalentContributionsAssebler>  equivalentContributionsAssemblers)
        {
            this.solver = solver;
            this.subdomainUpdaters = subdomainUpdaters;
            this.mappings = mappings;
            this.linearSystems = linearSystems;
            this.provider = provider;
            this.increments = increments;
            this.totalDOFs = totalDOFs;
            this.globalRHS = new Vector(totalDOFs);
            this.u = uInitialFreeDOFDisplacementsPerSubdomain; //prosthiki MS, TODO:possibly check for compatibility elements format: u.Add(subdomain.ID, new Vector(subdomain.RHS.Length));
            //this.uPlusdu = uInitialFreeDOFDisplacementsPerSubdomain; //prosthiki MS: commented out possible pass by reference
            this.boundaryNodes = boundaryNodes;
            this.initialConvergedBoundaryDisplacements = initialConvergedBoundaryDisplacements;
            this.totalBoundaryDisplacements = totalBoundaryDisplacements;
            this.equivalentContributionsAssemblers = equivalentContributionsAssemblers;

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
            //u.Clear(); prosthiki MS
            du.Clear();
            uPlusdu.Clear(); //prosthiki MS

            foreach (ILinearSystem subdomain in linearSystems)
            {
                Vector r = new Vector(subdomain.RHS.Length);
                ((Vector)subdomain.RHS).CopyTo(r.Data, 0);
                r.Multiply(1 / (double)increments); 
                rhs.Add(subdomain.ID, r); //comment MS: xtizei to rhs, field tou NRNLAnalyzer ok, swsta giati rhs sth dhmiourgia twn linear systems exoume perasei to model.Subdomains[0].Forces
                //u.Add(subdomain.ID, new Vector(subdomain.RHS.Length)); //prosthiki MS
                du.Add(subdomain.ID, new Vector(subdomain.RHS.Length));
                Vector tempCopy = new Vector(u[subdomain.ID].Length);
                u[subdomain.ID].Data.CopyTo(tempCopy.Data, 0);
                uPlusdu.Add(subdomain.ID, tempCopy); //prosthiki MS
                mappings[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].SubdomainToGlobalVector(((Vector)subdomain.RHS).Data, globalRHS.Data);
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
                rhs[subdomain.ID] = r; //comment MS: xanahtizei to rhs kai to xanaperna sto globalRHS(molis to kane clear opote ok)
                mappings[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].SubdomainToGlobalVector(((Vector)subdomain.RHS).Data, globalRHS.Data);
            }
            rhsNorm = provider.RHSNorm(globalRHS.Data);
        }

        public void Initialize()
        {
            solver.Initialize();
        }

        private void UpdateRHS(int step)
        {
            foreach (ILinearSystem subdomain in linearSystems)
            {
                Vector subdomainRHS = ((Vector)subdomain.RHS);
                rhs[subdomain.ID].CopyTo(subdomainRHS.Data, 0);
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
                UpdateRHS(increment);//comment MS2: apo to rhs[subdomain.ID] pernaei sto subdomain.RHS h fixed timh (externalLoads/increments) (ginetai copy kai oxi add)  AFTO thewreitai RHS sthn prwth iteration
                UpdateRHSForLinearizationContributions(increment + 1, increments);// opote edw pouu uparxei to neo rhs tou kuklou epanalhpsewn prepei na afairountai oi draseis
                double firstError = 0;
                int step = 0;
                for (step = 0; step < maxSteps; step++)
                {
                    solver.Solve();
                    errorNorm = rhsNorm != 0 ? CalculateInternalRHS(increment, step,increments) / rhsNorm : 0;//comment MS2: to subdomain.RHS lamvanei thn timh nIncrement*(externalLoads/increments)-interanalRHS me xrhsh ths fixed timhs apo to rhs[subdomain.ID]
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

            StoreLogResults(start, end);
        }

        private double CalculateInternalRHS(int currentIncrement, int step, int totalIncrements) 
        {
            globalRHS.Clear();
            foreach (ILinearSystem subdomain in linearSystems)
            {
                if (currentIncrement == 0 && step == 0)
                {
                    //Array.Clear(du[subdomain.ID].Data, 0, du[subdomain.ID].Length);
                    //Array.Clear(uPlusdu[subdomain.ID].Data, 0, uPlusdu[subdomain.ID].Length);
                    //du[subdomain.ID].Add(((Vector)subdomain.Solution));
                    //uPlusdu[subdomain.ID].Add(((Vector)subdomain.Solution));
                    //du[subdomain.ID].Subtract(u[subdomain.ID]);
                    du[subdomain.ID].Add(((Vector)subdomain.Solution));
                    Array.Clear(uPlusdu[subdomain.ID].Data, 0, uPlusdu[subdomain.ID].Length);
                    uPlusdu[subdomain.ID].Add(u[subdomain.ID]);
                    uPlusdu[subdomain.ID].Add(du[subdomain.ID]);
                }
                else
                {
                    du[subdomain.ID].Add(((Vector)subdomain.Solution));
                    Array.Clear(uPlusdu[subdomain.ID].Data, 0, uPlusdu[subdomain.ID].Length);
                    uPlusdu[subdomain.ID].Add(u[subdomain.ID]);
                    uPlusdu[subdomain.ID].Add(du[subdomain.ID]);
                }
                //Vector<double> internalRHS = (Vector<double>)subdomain.GetRHSFromSolution(u[subdomain.ID], du[subdomain.ID]);
                Vector internalRHS = (Vector)subdomainUpdaters[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].GetRHSFromSolutionWithInitialDisplacemntsEffect(uPlusdu[subdomain.ID], du[subdomain.ID],  boundaryNodes,
                initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, currentIncrement+1, totalIncrements);//TODOMaria this calculates the internal forces
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

        private void UpdateSolution()
        {
            foreach (ILinearSystem subdomain in linearSystems)
            {
                u[subdomain.ID].Add(du[subdomain.ID]);
            }
        }

        public void SaveMaterialState()
        {
            foreach (ILinearSystem subdomain in linearSystems)
            {
                subdomainUpdaters[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].UpdateState();
            }
        }

        public Dictionary<int, Vector> GetConvergedSolutionVectorsOfFreeDofs()
        {
            return u;
            // return uplusdu einai to idio afou exei ginei molis UpdateSolution
        }

        public Dictionary<int, Vector> GetConvergedIncrementalSolutionVectorsOfFreeDofs()
        {
            return du;
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

        private void UpdateRHSForLinearizationContributions(int nIncrement, int increments)
        {
            globalRHS.Clear();
            foreach (ILinearSystem subdomain in linearSystems)
            {
                //var equivalentContributionsAssembler = equivalentContributionsAssemblers[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index];
                var equivalentContributionsAssembler = equivalentContributionsAssemblers[subdomain.ID];
                Vector Contribution = equivalentContributionsAssembler.CalculateKfreeprescribedUpMultiplicationForSubdRHSContribution(
                    boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, increments);
                Vector subdomainRHS = ((Vector)subdomain.RHS);
                subdomainRHS.Subtract(Contribution);
                mappings[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].SubdomainToGlobalVector(subdomainRHS.Data, globalRHS.Data);
            }
            rhsNorm = provider.RHSNorm(globalRHS.Data);
            //TODOMS: possibly it is not nessesary to update globalRHS (and of course not clear it before the loop) as it is not updated in 177-UpdateRHS(increment) either.
            // and it is used when it is recalculated in CalculateInternalRHS......
        }

        #endregion
    }
}
