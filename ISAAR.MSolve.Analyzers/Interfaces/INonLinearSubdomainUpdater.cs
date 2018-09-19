using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System.Collections.Generic;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface INonLinearSubdomainUpdater
    {
        void ScaleConstraints(double scalingFactor);
        IVector GetRHSFromSolution(IVector solution, IVector dSolution);
        void UpdateState();
        void ResetState();
    }
}
