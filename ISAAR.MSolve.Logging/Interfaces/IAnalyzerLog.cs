using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.Logging.Interfaces
{
    public interface IAnalyzerLog
    {
        void StoreResults(DateTime startTime, DateTime endTime, IVector<double> solution);
    }
}
