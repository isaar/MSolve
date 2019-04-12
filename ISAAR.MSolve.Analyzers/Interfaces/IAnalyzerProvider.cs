namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IAnalyzerProvider
    {
        //TODO: This should be accessed by the solver. Any element matrix providers should be passed there.
        IDirichletEquivalentLoadsAssembler DirichletLoadsAssembler { get; }

        void Reset();
    }
}
