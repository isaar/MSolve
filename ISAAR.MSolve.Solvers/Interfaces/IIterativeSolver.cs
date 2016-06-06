using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface IIterativeSolver : ISolver
    {
        int CurrentIteration { get; }
        void Initialize(IVector<double> x, IVector<double> residual, double detf);
        void Solve(int maxIterations, double tolerance);
    }
}
