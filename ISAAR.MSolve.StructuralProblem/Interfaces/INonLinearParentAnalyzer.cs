using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface INonLinearParentAnalyzer : IAnalyzer
    {
        double[] GetOtherRHSComponents(ISolverSubdomain subdomain, IVector<double> currentSolution);
    }
}
