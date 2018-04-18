namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface INonLinearProvider : IAnalyzerProvider
    {
        // TODO: ID fields should be removed and the whole domain decomposition thing must be decoupled from the analyzers
        double RHSNorm(double[] rhs);
        void ProcessInternalRHS(int id, double[] rhs, double[] solution);
    }
}
