using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Analyzers.NonLinear
{
    //TODO: Confusing name. The child analyzer of this is a nonlinear analyzer.
    public interface INonLinearParentAnalyzer_v2 : IParentAnalyzer
    {
        IVector GetOtherRhsComponents(ILinearSystem_v2 linearSystem, IVector currentSolution);
    }
}
