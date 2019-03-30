using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Analyzers
{
    public class LinearAnalyzer_v2 : IChildAnalyzer
    {
        private readonly IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems;
        private readonly IStructuralModel_v2 model;
        private readonly IAnalyzerProvider_v2 provider;
        private readonly ISolver_v2 solver;

        public LinearAnalyzer_v2(IStructuralModel_v2 model, ISolver_v2 solver, IAnalyzerProvider_v2 provider)
        {
            this.model = model;
            this.solver = solver;
            this.linearSystems = solver.LinearSystems;
            this.provider = provider;
        }

        public Dictionary<int, ILogFactory_v2> LogFactories { get; } = new Dictionary<int, ILogFactory_v2>();
        public Dictionary<int, IAnalyzerLog_v2[]> Logs { get; } = new Dictionary<int, IAnalyzerLog_v2[]>();

        public IParentAnalyzer ParentAnalyzer { get; set; }

        public void BuildMatrices()
        {
            if (ParentAnalyzer == null) throw new InvalidOperationException("This linear analyzer has no parent.");

            ParentAnalyzer.BuildMatrices();
            //solver.Initialize();
        }

        public void Initialize(bool isFirstAnalysis)
        {
            InitializeLogs();
            solver.Initialize();
        }

        public void Solve()
        {
            DateTime start = DateTime.Now;
            AddEquivalentNodalLoadsToRHS(); //TODO: The initial rhs (from other loads) should also be built by the analyzer instead of the model.
            solver.Solve();
            DateTime end = DateTime.Now;
            StoreLogResults(start, end);
        }

        private void AddEquivalentNodalLoadsToRHS()
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                try
                {
                    // Make sure there is at least one non zero prescribed displacement.
                    (INode node, DOFType dof, double displacement) = linearSystem.Subdomain.Constraints.Find(du => du != 0.0);

                    //TODO: the following 2 lines are meaningless outside diplacement control (and even then, they are not so clear).
                    double scalingFactor = 1;
                    IVector initialFreeSolution = linearSystem.CreateZeroVector();

                    IVector equivalentNodalLoads = provider.DirichletLoadsAssembler.GetEquivalentNodalLoads(
                        linearSystem.Subdomain, initialFreeSolution, scalingFactor);
                    linearSystem.RhsVector.SubtractIntoThis(equivalentNodalLoads);
                }
                catch (KeyNotFoundException)
                {
                    // There aren't any non zero prescribed displacements, therefore we do not have to calculate the equivalent 
                    // nodal loads, which is an expensive operation (all elements are accessed, their stiffness is built, etc..)
                }
            }
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
                    l.StoreResults(start, end, linearSystems[id].Solution);
        }
    }
}
