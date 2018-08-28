using System;
using System.Collections.Generic;
using System.Text;

//TODO: The logs should be created by the user and injected into the analyzers. This way they are accessible by the
//      user during and after the analysis.
namespace ISAAR.MSolve.Logging.Interfaces
{
    /// <summary>
    /// Used by the analyzers to create the logs.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface ILogFactory
    {
        IAnalyzerLog[] CreateLogs();
    }
}
