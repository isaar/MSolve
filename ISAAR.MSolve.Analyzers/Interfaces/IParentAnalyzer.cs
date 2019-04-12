namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IParentAnalyzer: IAnalyzer
    {
        /// <summary>
        /// The child analyzer should be injected in the constructor, since it is mandatory for the parent analyzer to work.
        /// </summary>
        IChildAnalyzer ChildAnalyzer { get; }
    }
}
