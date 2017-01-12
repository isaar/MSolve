using System.Linq;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.PCG
{
    public class SolverPCGDiagonalPreconditioner : IIterativeSolverPreconditioner
    {
        #region ISolverPCGMatrixCalculator Members

        private readonly ILinearSystem linearSystem;
        private readonly double[] diagonalPreconditioner;

        public SolverPCGDiagonalPreconditioner(ILinearSystem linearSystem)
        {
            this.linearSystem = linearSystem;
            diagonalPreconditioner = new double[linearSystem.Matrix.Rows];
            for (int i = 0; i < diagonalPreconditioner.Length; i++)
                diagonalPreconditioner[i] = 1d / linearSystem.Matrix[i, i];
        }

        public void Precondition(IVector vIn, IVector vOut)
        {
            for (int i = 0; i < vIn.Length; i++)
                vOut[i] = diagonalPreconditioner[i] * vIn[i];
        }

        #endregion
    }
}
