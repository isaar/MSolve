namespace ISAAR.MSolve.Analyzers
{
    public interface IParentAnalyzer: IAnalyzer_v2
    {
        /// <summary>
        /// The child analyzer should be injected in the constructor, since it is mandatory for the parent analyzer to work.
        /// </summary>
        IChildAnalyzer ChildAnalyzer { get; }
    }
}
