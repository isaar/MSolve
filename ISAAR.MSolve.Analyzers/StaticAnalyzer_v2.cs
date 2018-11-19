using System;
using System.Collections.Generic;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Analyzers
{
    public class StaticAnalyzer_v2 : IAnalyzer_v2, INonLinearParentAnalyzer_v2
    {
        private readonly IReadOnlyList<ILinearSystem_v2> linearSystems;
        private readonly ISolver_v2 solver;
        private readonly IStaticProvider_v2 provider;
        private IAnalyzer_v2 childAnalyzer;
        private IAnalyzer_v2 parentAnalyzer = null;
        private readonly Dictionary<int, IAnalyzerLog[]> logs = new Dictionary<int, IAnalyzerLog[]>();
        private bool areDofsOrdered = false;

        public StaticAnalyzer_v2(ISolver_v2 solver, IStaticProvider_v2 provider, IAnalyzer_v2 embeddedAnalyzer)
        {
            this.solver = solver;
            this.linearSystems = solver.LinearSystems;
            this.provider = provider;
            this.childAnalyzer = embeddedAnalyzer;
            this.childAnalyzer.ParentAnalyzer = this;
        }

        private void InitalizeMatrices()
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems) provider.CalculateMatrix(linearSystem);

            //provider.CalculateMatrices();
            //subdomain.Matrix = provider.Ks[subdomain.ID];
        }

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog[]> Logs { get { return logs; } }

        public IAnalyzer_v2 ParentAnalyzer
        {
            get { return parentAnalyzer; }
            set { parentAnalyzer = value; }
        }

        public IAnalyzer_v2 ChildAnalyzer
        {
            get { return childAnalyzer; }
            set { childAnalyzer = value; }
        }

        public void Initialize()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Static analyzer must contain an embedded analyzer.");
            //InitalizeMatrices();
            childAnalyzer.Initialize();
        }

        public void Solve()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Static analyzer must contain an embedded analyzer.");
            childAnalyzer.Solve();
        }

        public void BuildMatrices()
        {
            //TODO: Dof ordering should be handled separately from matrix building. It should be done in Initialize() (for this
            //      analyzer), which should be called before BuildMatrices(). Actually BuildMatrices() should not be called by
            //      the user.
            if (!areDofsOrdered)
            {
                solver.OrderDofs();
                areDofsOrdered = true;
            }

            InitalizeMatrices();
            //int rows = subdomains[1].Matrix.Rows;
            //double[,] data = new double[rows, rows];
            //for (int i = 0; i < rows; i++)
            //    for (int j = 0; j < rows; j++)
            //        data[i, j] = subdomains[1].Matrix[i, j];

            //var m = new Matrix2D<double>(data);
            //var w = new double[rows];
            //var v = new double[rows, rows];
            //m.SVD(w, v);
        }

        #endregion

        #region INonLinearParentAnalyzer Members

        public IVector GetOtherRHSComponents(ILinearSystem_v2 linearSystem, IVector currentSolution)
        {
            //TODO: use a ZeroVector class that avoid doing useless operations or refactor this method. E.g. let this method 
            // alter the child analyzer's rhs vector, instead of the opposite (which is currently done).
            return linearSystem.CreateZeroVector();
        }

        #endregion
    }
}
