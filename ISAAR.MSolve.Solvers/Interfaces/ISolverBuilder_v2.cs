using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface ISolverBuilder_v2
    {
        IDofOrderer DofOrderer { get; set; }
        ISolver_v2 BuildSolver(IStructuralModel_v2 model);
        ISolverBuilder_v2 Clone();
    }
}
