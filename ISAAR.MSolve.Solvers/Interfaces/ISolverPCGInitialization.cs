using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface ISolverPCGInitialization
    {
        double InitializeAndGetResidual(IList<ISolverSubdomain> subdomains, double[] r, double[] x);
    }
}
