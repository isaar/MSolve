using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface INonLinearParentAnalyzer : IAnalyzer
    {
        // TODO: ID fields should be removed and the whole domain decomposition thing must be decoupled from the analyzers
        double[] GetOtherRHSComponents(int id, IVector currentSolution);
    }
}
