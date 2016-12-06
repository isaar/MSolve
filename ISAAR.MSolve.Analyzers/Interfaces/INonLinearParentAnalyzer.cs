using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface INonLinearParentAnalyzer : IAnalyzer
    {
        double[] GetOtherRHSComponents(ILinearSystem subdomain, IVector currentSolution);
    }
}
