using ISAAR.MSolve.Logging.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

//TODO: should child analyzers hold references to their parent analyzers?
namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IAnalyzer_v2
    {
        Dictionary<int, IAnalyzerLog[]> Logs { get; }

        void BuildMatrices(); //This makes sense for parent analyzers only.
        void Initialize(); // The user should not have to call this.
        void Solve();
    }
}
