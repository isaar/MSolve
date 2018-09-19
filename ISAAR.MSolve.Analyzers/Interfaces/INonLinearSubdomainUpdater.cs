using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System.Collections.Generic;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface INonLinearSubdomainUpdater
    {
        void ApplyConstraints(Dictionary<int, Dictionary<DOFType, double>> globalConstraintsDictionary);
        IVector GetRHSFromSolution(IVector solution, IVector dSolution);
        void UpdateState();
        void ResetState();
    }
}
