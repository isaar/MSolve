using System.Collections.Generic;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface IMatrixSolver
    {
        IDictionary<int, IMatrixLinearSystem> MatrixLinearSystems { get; }
    }
}
