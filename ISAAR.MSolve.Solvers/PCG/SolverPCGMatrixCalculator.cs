using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.PCG
{
    public class SolverPCGMatrixCalculator : IIterativeSolverMatrixCalculator
    {
        #region ISolverPCGMatrixCalculator Members

        private readonly IMatrixLinearSystem linearSystem;

        public SolverPCGMatrixCalculator(IMatrixLinearSystem linearSystem)
        {
            this.linearSystem = linearSystem;
        }

        public void MultiplyWithMatrix(IVectorOLD vIn, IVectorOLD vOut)
        {
            linearSystem.Matrix.Multiply(vIn, ((Vector)vOut).Data);
        }

        #endregion
    }
}
