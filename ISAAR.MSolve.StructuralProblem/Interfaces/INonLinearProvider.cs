using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface INonLinearProvider : IAnalyzerProvider
    {
        double RHSNorm(double[] rhs);
        void ProcessInternalRHS(ISolverSubdomain subdomain, double[] rhs, double[] solution);
    }
}
