using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System.Collections.Generic;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface IIterativeSolverInitializer
    {
        double InitializeAndGetResidual(IVectorOLD r, IVectorOLD x);
    }
}
