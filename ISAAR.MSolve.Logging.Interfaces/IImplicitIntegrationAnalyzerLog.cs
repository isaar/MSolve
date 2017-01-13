using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Logging.Interfaces
{
    public interface IImplicitIntegrationAnalyzerLog
    {
        void StoreResults(DateTime startTime, DateTime endTime, IAnalyzerLog log);
    }
}
