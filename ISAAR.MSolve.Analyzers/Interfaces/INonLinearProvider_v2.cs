using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface INonLinearProvider_v2 : IAnalyzerProvider_v2
    {
        double CalculateRhsNorm(IVectorView rhs);
        void ProcessInternalRhs(ILinearSystem_v2 linearSystem, IVector rhs, IVectorView solution); //TODO: this does nothing
    }
}
