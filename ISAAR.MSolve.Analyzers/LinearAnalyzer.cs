using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    public class LinearAnalyzer : IAnalyzer
    {
        private IAnalyzer parentAnalyzer = null;
        private readonly ISolver solver;
        private readonly IDictionary<int, ILinearSystem> linearSystems;
        private readonly Dictionary<int, ILogFactory> logFactories = new Dictionary<int, ILogFactory>();
        private readonly Dictionary<int, IAnalyzerLog[]> logs = new Dictionary<int, IAnalyzerLog[]>();

        public LinearAnalyzer(ISolver solver, IDictionary<int, ILinearSystem> linearSystems)
        {
            this.solver = solver;
            this.linearSystems = linearSystems;
        }

        //TODO: this is temporarily injected as a property, so that the users only need to provide it if they need it.
        //      Normally it should be defined by the model depending on the boundary & initial conditions, body forces, etc.
        //      The provider (e.g. problem structural) may also need to be taken into consideration.
        public IDictionary<int, IEquivalentLoadsAssembler> EquivalentLoadsAssemblers { get; set; }

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

        public Dictionary<int, ILogFactory> LogFactories { get { return logFactories; } }
        public ISolver Solver { get { return solver; } }

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog[]> Logs { get { return logs; } }
        public IAnalyzer ParentAnalyzer
        {
            get { return parentAnalyzer; }
            set { parentAnalyzer = value; }
        }

        public IAnalyzer ChildAnalyzer
        {
            get { return null; }
            set {  throw new InvalidOperationException("Linear analyzer cannot contain an embedded analyzer."); }
        }

        public void Initialize()
        {
            InitializeLogs();
            solver.Initialize();
        }

        public void Solve()
        {
            DateTime start = DateTime.Now;
            AddEquivalentNodalLoadsToRHS(); //TODO: The initial rhs should also be built by the analyzer instead of the model.
            solver.Solve();
            DateTime end = DateTime.Now;
            StoreLogResults(start, end); 
        }

        public void BuildMatrices()
        {
            if (parentAnalyzer == null) throw new InvalidOperationException("This linear analyzer has no parent.");

            parentAnalyzer.BuildMatrices();
            //solver.Initialize();
        }

        #endregion

        private void AddEquivalentNodalLoadsToRHS()
        {
            if (EquivalentLoadsAssemblers != null)
            {
                foreach (var subdomain in linearSystems.Keys)
                {
                    Vector subdomainRHS = ((Vector)linearSystems[subdomain].RHS);
                    var equivalentLoadsAssembler = EquivalentLoadsAssemblers[subdomain];

                    var initialFreeSolution = new Vector(subdomainRHS.Length);
                    double scalingFactor = 1;
                    var equivalentNodalLoads = (Vector)equivalentLoadsAssembler.GetEquivalentNodalLoads(
                        initialFreeSolution, scalingFactor);
                    subdomainRHS.Subtract(equivalentNodalLoads);
                }
            }
        }
    }
}
