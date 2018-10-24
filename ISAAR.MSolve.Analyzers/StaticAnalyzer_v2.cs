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
        private readonly IDictionary<int, ILinearSystem_v2> subdomains;
        private readonly IStaticProvider_v2 provider;
        private IAnalyzer_v2 childAnalyzer;
        private IAnalyzer_v2 parentAnalyzer = null;
        private readonly Dictionary<int, IAnalyzerLog[]> logs = new Dictionary<int, IAnalyzerLog[]>();

        public StaticAnalyzer_v2(IStaticProvider_v2 provider, IAnalyzer_v2 embeddedAnalyzer, IDictionary<int, ILinearSystem_v2> subdomains)
        {
            this.provider = provider;
            this.childAnalyzer = embeddedAnalyzer;
            this.subdomains = subdomains;
            this.childAnalyzer.ParentAnalyzer = this;
        }

        private void InitalizeMatrices()
        {
            foreach (ILinearSystem_v2 subdomain in subdomains.Values)
                provider.CalculateMatrix(subdomain);
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

        public IVector GetOtherRHSComponents(ILinearSystem_v2 subdomain, IVector currentSolution)
        {
            return Vector.CreateZero(subdomain.RhsVector.Length);
        }

        #endregion
    }
}
