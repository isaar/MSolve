using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Analyzers
{
    public enum CrackPropagationTermination
    {
        RequiredIterationsWereCompleted, CrackExitsDomainBoundary, MechanismIsCreated, FractureToughnessIsExceeded
    }
}
