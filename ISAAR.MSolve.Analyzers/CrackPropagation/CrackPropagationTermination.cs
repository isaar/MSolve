using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Analyzers.CrackPropagation
{
    public enum CrackPropagationTermination
    {
        RequiredIterationsWereCompleted, CrackExitsDomainBoundary, MechanismIsCreated, FractureToughnessIsExceeded
    }
}
