using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    //TODO: Confusing name. The child analyzer of this is a nonlinear analyzer.
    public interface INonLinearParentAnalyzer : IParentAnalyzer
    {
        IVector GetOtherRhsComponents(ILinearSystem linearSystem, IVector currentSolution);
    }
}
