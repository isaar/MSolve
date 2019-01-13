namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface ISolver
    {
        // TODO: I think this needs to be called only once (by the analyzer), not every time the matrix must be factorized.
        void Initialize();
        void Solve();
    }
}
