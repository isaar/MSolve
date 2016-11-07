using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface ISolverPCGMatrixCalculator
    {
        int VectorSize { get; }
        void Precondition(IVector vIn, IVector vOut);
        void MultiplyWithMatrix(IVector vIn, IVector vOut);
    }
}
