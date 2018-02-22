namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface ISubdomainGlobalMapping
    {
        void SubdomainToGlobalVector(double[] vIn, double[] vOut);
        void SubdomainToGlobalVectorMeanValue(double[] vIn, double[] vOut);
        void SplitGlobalVectorToSubdomain(double[] vIn, double[] vOut);
    }
}
