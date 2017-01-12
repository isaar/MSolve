using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface INonLinearSubdomainUpdater
    {
        IVector GetRHSFromSolution(IVector solution, IVector dSolution);
        void UpdateState();
        void ResetState();
    }
}
