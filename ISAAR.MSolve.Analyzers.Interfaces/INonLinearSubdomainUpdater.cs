using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface INonLinearSubdomainUpdater
    {
        IVectorOLD GetRHSFromSolution(IVectorOLD solution, IVectorOLD dSolution);
        void UpdateState();
        void ResetState();
    }
}
