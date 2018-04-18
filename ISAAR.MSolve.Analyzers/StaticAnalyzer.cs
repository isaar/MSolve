using System;
using System.Collections.Generic;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    public class StaticAnalyzer : IAnalyzer, INonLinearParentAnalyzer
    {
        private readonly IDictionary<int, IMatrixLinearSystem> subdomains;
        private readonly IStaticProvider provider;
        private IAnalyzer childAnalyzer;
        private IAnalyzer parentAnalyzer = null;
        private readonly Dictionary<int, IAnalyzerLog[]> logs = new Dictionary<int, IAnalyzerLog[]>();

        public StaticAnalyzer(IStaticProvider provider, IAnalyzer embeddedAnalyzer, IMatrixLinearSystem subdomain)
        {
            this.provider = provider;
            this.childAnalyzer = embeddedAnalyzer;
            this.childAnalyzer.ParentAnalyzer = this;
            this.subdomains = new Dictionary<int, IMatrixLinearSystem>();
            this.subdomains.Add(subdomain.ID, subdomain);
        }

        public StaticAnalyzer(IStaticProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, IMatrixLinearSystem> subdomains)
        {
            this.provider = provider;
            this.childAnalyzer = embeddedAnalyzer;
            this.subdomains = subdomains;
            this.childAnalyzer.ParentAnalyzer = this;
        }

        private void InitalizeMatrices()
        {
            foreach (IMatrixLinearSystem subdomain in subdomains.Values)
                subdomain.Matrix = provider.CalculateMatrix(subdomain.ID);
            //provider.CalculateMatrices();
                //subdomain.Matrix = provider.Ks[subdomain.ID];
        }

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog[]> Logs { get { return logs; } }
        public IAnalyzer ParentAnalyzer
        {
            get { return parentAnalyzer; }
            set { parentAnalyzer = value; }
        }

        public IAnalyzer ChildAnalyzer
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

        public double[] GetOtherRHSComponents(int id, IVector currentSolution)
        {
            return new double[currentSolution.Length];
        }

        #endregion
    }
}
