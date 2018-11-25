using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface INonLinearParentAnalyzer_v2 : IParentAnalyzer
    {
        IVector GetOtherRhsComponents(ILinearSystem_v2 linearSystem, IVector currentSolution);
    }
}
