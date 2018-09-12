using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Solvers.PCG
{
    public interface IPreconditionerBuilder
    {
        IPreconditioner BuildPreconditioner(IMatrixView matrix);
    }
}
