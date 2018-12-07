using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Analyzers
{
    public class LinearAnalyzer_v2 : IChildAnalyzer
    {
        private readonly IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems;
        private readonly ISolver_v2 solver;

        public LinearAnalyzer_v2(ISolver_v2 solver)
        {
            this.solver = solver;
            this.linearSystems = solver.LinearSystems;
        }

        public Dictionary<int, ILogFactory> LogFactories { get; } = new Dictionary<int, ILogFactory>();
        public Dictionary<int, IAnalyzerLog[]> Logs { get; } = new Dictionary<int, IAnalyzerLog[]>();

        public IParentAnalyzer ParentAnalyzer { get; set; }

        public void BuildMatrices()
        {
            if (ParentAnalyzer == null) throw new InvalidOperationException("This linear analyzer has no parent.");

            ParentAnalyzer.BuildMatrices();
            //solver.Initialize();
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

        private void InitializeLogs()
        {
            Logs.Clear();
            foreach (int id in LogFactories.Keys) Logs.Add(id, LogFactories[id].CreateLogs());
        }

        private void StoreLogResults(DateTime start, DateTime end)
        {
            foreach (int id in Logs.Keys)
                foreach (var l in Logs[id])
                    l.StoreResults(start, end, linearSystems[id].Solution.ToLegacyVector());
        }
    }
}
