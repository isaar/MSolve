using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.Solvers.PCG
{
    public class SolverPCGSimpleSearchVectorCalculator : IPCGSearchVectorCalculator
    {
        private int searchVectorSize = 0;
        private Vector<double> zOld, rOld;

        #region ISearchVectorCalculator Members

        public void CalculateSearchVector(IIterativeSolver solver)
        {
            SolverPCG pcg = (SolverPCG)solver;
            if (pcg.CurrentIteration > 0)
            {
                double b = (pcg.VectorZ * pcg.VectorR) / (zOld * rOld);
                for (int i = 0; i < searchVectorSize; i++) pcg.VectorP[i] = pcg.VectorZ[i] + b * pcg.VectorP[i];
            }
            else
            {
                searchVectorSize = pcg.VectorZ.Length;
                zOld = new Vector<double>(searchVectorSize);
                rOld = new Vector<double>(searchVectorSize);
                pcg.VectorZ.CopyTo(pcg.VectorP.Data, 0);
            }
            pcg.VectorZ.CopyTo(zOld.Data, 0);
            pcg.VectorR.CopyTo(rOld.Data, 0);
        }

        public double CalculateGradient(IIterativeSolver solver)
        {
            SolverPCG pcg = (SolverPCG)solver;
            return (pcg.VectorZ * pcg.VectorR) / (pcg.VectorP * pcg.VectorQ);
        }

        public bool InitializeStartingVectorFromSearchVectors(IVector<double> x, IVector<double> b)
        {
            return false;
        }

        public void ClearSearchVectors(int vectorsToKeepFromTop)
        {
        }

        #endregion
    }
}
