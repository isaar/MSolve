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

        //TODO: The following two properties make sense for either parent or child analyzers. Perhaps a similar effect could be 
        //      achieved cleanly with decorators. For now use 2 interfaces: IParentAnalyzer and IChildAnalyzer
        IAnalyzer_v2 ChildAnalyzer { get; set; }
        IAnalyzer_v2 ParentAnalyzer { get; set; }

        void BuildMatrices(); //This makes sense for parent analyzers only.
        void Initialize(); // The user should not have to call this.
        void Solve();
    }
}
