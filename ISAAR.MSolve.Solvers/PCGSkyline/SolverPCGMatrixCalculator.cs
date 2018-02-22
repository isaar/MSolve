using System.Linq;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.PCGSkyline
{
    public class SolverPCGMatrixCalculator<T> : ISolverPCGMatrixCalculator where T : IMatrix2D
    {
        #region ISolverPCGMatrixCalculator Members

        private readonly SolverPCG<T> solver;

        public SolverPCGMatrixCalculator(SolverPCG<T> solver)
        {
            this.solver = solver;
        }

        public int VectorSize 
        { 
            get 
            { 
                return solver.LinearSystem.RHS.Length; 
            } 
        }

        public void Precondition(IVector vIn, IVector vOut)
        {
            solver.Precondition(vIn, vOut);
        }

        public void MultiplyWithMatrix(IVector vIn, IVector vOut)
        {
            solver.MultiplyWithMatrix(vIn, vOut);
        }

        #endregion
    }
}
