using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface ILinearSubdomainUpdater
    {
        IVector GetRHS();
    }
}
