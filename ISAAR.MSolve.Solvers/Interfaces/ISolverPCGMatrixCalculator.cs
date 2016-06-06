using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface ISolverPCGMatrixCalculator
    {
        int VectorSize { get; }
        void Precondition(IVector<double> vIn, IVector<double> vOut);
        void MultiplyWithMatrix(IVector<double> vIn, IVector<double> vOut);
    }
}
