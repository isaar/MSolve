using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.Solvers.PCG
{
    public class SolverPCGReorthogonalizationSearchVectorCalculator : IPCGSearchVectorCalculator
    {
        private readonly List<Vector<double>> p = new List<Vector<double>>();
        private readonly List<Vector<double>> q = new List<Vector<double>>();

        #region ISearchVectorCalculator Members

        public void CalculateSearchVector(IIterativeSolver solver)
        {
            SolverPCG pcg = (SolverPCG)solver;
            SearchVectors.CalculateReorthogonalizedSearchVector(pcg.VectorZ, pcg.VectorP, p, q);
        }

        public double CalculateGradient(IIterativeSolver solver)
        {
            SolverPCG pcg = (SolverPCG)solver;
            return SearchVectors.CalculateReorthogonalizedGradient(pcg.VectorP, pcg.VectorQ, pcg.VectorR, p, q);
        }

        public bool InitializeStartingVectorFromSearchVectors(IVector<double> x, IVector<double> b)
        {
            return SearchVectors.InitializeStartingVectorFromReorthogonalizedSearchVectors(
                (Vector<double>)x, (Vector<double>)b, p, q);
        }

        public void ClearSearchVectors(int vectorsToKeepFromTop)
        {
            if (vectorsToKeepFromTop < 1)
            {
                p.Clear();
                q.Clear();
            }
            else
            {
                var newps = new Vector<double>[vectorsToKeepFromTop];
                var newqs = new Vector<double>[vectorsToKeepFromTop];
                for (int i = 0; i < vectorsToKeepFromTop; i++)
                {
                    newps[i] = p[i];
                    newqs[i] = q[i];
                }

                p.Clear();
                q.Clear();

                for (int i = 0; i < vectorsToKeepFromTop; i++)
                {
                    p.Add(newps[i]);
                    q.Add(newqs[i]);
                }
            }
        }

        #endregion
    }
}
