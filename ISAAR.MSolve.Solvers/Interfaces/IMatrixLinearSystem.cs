using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface IMatrixLinearSystem
    {
        int ID { get; }
        IMatrix2D Matrix { get; set; }
        IVector RHS { get; }
        IVector Solution { get; set;  }
    }
}
