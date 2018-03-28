using System.Collections.Generic;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.Solvers.PCG
{
    public class SolverPCGReorthogonalizationSearchVectorCalculator : IPCGSearchVectorCalculator
    {
        private readonly List<IVectorOLD> p = new List<IVectorOLD>();
        private readonly List<IVectorOLD> q = new List<IVectorOLD>();

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

        public bool InitializeStartingVectorFromSearchVectors(IVectorOLD x, IVectorOLD b)
        {
            return SearchVectors.InitializeStartingVectorFromReorthogonalizedSearchVectors(x, b, p, q);
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
                var newps = new IVectorOLD[vectorsToKeepFromTop];
                var newqs = new IVectorOLD[vectorsToKeepFromTop];
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
