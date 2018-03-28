using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface IMatrixLinearSystem
    {
        int ID { get; }
        IMatrix2D Matrix { get; set; }
        IVectorOLD RHS { get; }
        IVectorOLD Solution { get; set;  }
    }
}
