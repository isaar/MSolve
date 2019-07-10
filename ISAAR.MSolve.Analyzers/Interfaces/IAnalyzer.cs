using System.Collections.Generic;
using ISAAR.MSolve.Logging.Interfaces;

//TODO: should child analyzers hold references to their parent analyzers?
namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IAnalyzer
    {
        Dictionary<int, IAnalyzerLog[]> Logs { get; }

        void BuildMatrices(); //This makes sense for parent analyzers only.
        void Initialize(bool isFirstAnalysis); // The user should not have to call this.
        void Solve();
    }
}
