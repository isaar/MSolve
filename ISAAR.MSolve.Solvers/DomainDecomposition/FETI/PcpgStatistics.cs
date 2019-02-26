using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.FETI
{
    internal class PcpgStatistics
    {
        internal bool HasConverged { get; set; }
        internal int NumIterations { get; set; }
        internal double ResidualNormEstimateRatio { get; set; }
    }
}
