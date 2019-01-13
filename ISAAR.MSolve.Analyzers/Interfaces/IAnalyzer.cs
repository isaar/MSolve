using ISAAR.MSolve.Logging.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IAnalyzer
    {
        IAnalyzer ParentAnalyzer { get; set; } // This makes sense for child analyzers only. Could it be done with decorators cleanly?
        IAnalyzer ChildAnalyzer { get; set; } // Ditto, but this makes sense for parent analyzers only.
        Dictionary<int, IAnalyzerLog[]> Logs { get; }
        void BuildMatrices(); //This makes sense for parent analyzers only. The user should not have to call this.
        void Initialize(); // The user should not have to call this.
        void Solve();
    }
}
