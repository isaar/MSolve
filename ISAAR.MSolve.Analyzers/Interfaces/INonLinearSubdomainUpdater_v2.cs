
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface INonLinearSubdomainUpdater_v2
    {
        void ScaleConstraints(double scalingFactor);
        IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution);
        void UpdateState();
        void ResetState();
    }
}
