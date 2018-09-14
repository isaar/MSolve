using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    public class LinearAnalyzer : IAnalyzer
    {
        private IAnalyzer parentAnalyzer = null;
        private readonly ISolver solver;
        private readonly IDictionary<int, ILinearSystem> subdomains;
        private readonly Dictionary<int, ILogFactory> logFactories = new Dictionary<int, ILogFactory>();
        private readonly Dictionary<int, IAnalyzerLog[]> logs = new Dictionary<int, IAnalyzerLog[]>();

        public LinearAnalyzer(ISolver solver, IDictionary<int, ILinearSystem> subdomains)
        {
            this.solver = solver;
            this.subdomains = subdomains;
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
                    l.StoreResults(start, end, subdomains[id].Solution);
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
    }
}
