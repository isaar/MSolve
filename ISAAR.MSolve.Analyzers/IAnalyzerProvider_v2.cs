namespace ISAAR.MSolve.Analyzers
{
    public interface IAnalyzerProvider_v2
    {
        //TODO: This should be accessed by the solver. Any element matrix providers should be passed there.
        IDirichletEquivalentLoadsAssembler DirichletLoadsAssembler { get; }

        void Reset();
    }
}
