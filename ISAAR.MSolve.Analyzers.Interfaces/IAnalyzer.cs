namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IAnalyzer
    {
        IAnalyzer ParentAnalyzer { get; set; }
        IAnalyzer ChildAnalyzer { get; set; }
        //Dictionary<int, IAnalyzerLog[]> Logs { get; }
        void BuildMatrices();
        void Initialize();
        void Solve();
    }
}
