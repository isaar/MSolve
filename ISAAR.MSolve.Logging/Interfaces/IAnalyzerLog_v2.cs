using ISAAR.MSolve.LinearAlgebra.Vectors;
using System;

namespace ISAAR.MSolve.Logging.Interfaces
{
    public interface IAnalyzerLog_v2
    {
        void StoreResults(DateTime startTime, DateTime endTime, IVectorView solution);
    }
}
