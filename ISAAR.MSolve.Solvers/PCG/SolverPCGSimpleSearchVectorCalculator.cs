using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;


namespace ISAAR.MSolve.Solvers.PCG
{
    public class SolverPCGSimpleSearchVectorCalculator : IPCGSearchVectorCalculator
    {
        private int searchVectorSize = 0;
        private IVector zOld, rOld;

        #region ISearchVectorCalculator Members

        public void CalculateSearchVector(IIterativeSolver solver)
        {
            SolverPCG pcg = (SolverPCG)solver;
            if (pcg.CurrentIteration > 0)
            {
                double b = pcg.VectorZ.DotProduct(pcg.VectorR) / zOld.DotProduct(rOld);
                for (int i = 0; i < searchVectorSize; i++) pcg.VectorP[i] = pcg.VectorZ[i] + b * pcg.VectorP[i];
            }
            else
            {
                searchVectorSize = pcg.VectorZ.Length;
                zOld = new Vector(searchVectorSize);
                rOld = new Vector(searchVectorSize);
                pcg.VectorP.CopyFrom(0, searchVectorSize, pcg.VectorZ, 0);
                //pcg.VectorZ.CopyTo(pcg.VectorP.Data, 0);
            }
            zOld.CopyFrom(0, searchVectorSize, pcg.VectorZ, 0);
            rOld.CopyFrom(0, searchVectorSize, pcg.VectorR, 0);
            //pcg.VectorZ.CopyTo(zOld.Data, 0);
            //pcg.VectorR.CopyTo(rOld.Data, 0);
        }

        public double CalculateGradient(IIterativeSolver solver)
        {
            SolverPCG pcg = (SolverPCG)solver;
            return pcg.VectorZ.DotProduct(pcg.VectorR) / pcg.VectorP.DotProduct(pcg.VectorQ);
            //return (pcg.VectorZ * pcg.VectorR) / (pcg.VectorP * pcg.VectorQ);
        }

        public bool InitializeStartingVectorFromSearchVectors(IVector x, IVector b)
        {
            return false;
        }

        public void ClearSearchVectors(int vectorsToKeepFromTop)
        {
        }

        #endregion
    }
}
