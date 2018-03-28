using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;

namespace ISAAR.MSolve.Logging.Interfaces
{
    public interface IAnalyzerLog
    {
        void StoreResults(DateTime startTime, DateTime endTime, IVectorOLD solution);
    }
}
