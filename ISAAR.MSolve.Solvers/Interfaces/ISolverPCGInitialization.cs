using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System.Collections.Generic;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface ISolverPCGInitialization
    {
        double InitializeAndGetResidual(IList<ILinearSystem> linearSystems, IVector r, IVector x);
    }
}
