using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    public class FetiLogger
    {
        public int NumExpandedDomainFreeDofs { get; set; }
        public int NumLagrangeMultipliers { get; set; }
        public int NumUniqueGlobalFreeDofs { get; set; }
        public int PcpgIterations { get; set; }
        public double PcpgResidualNormEstimateRatio { get; set; }
    }
}
