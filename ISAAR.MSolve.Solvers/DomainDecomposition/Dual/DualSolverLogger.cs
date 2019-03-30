using System;
using System.Collections.Generic;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual
{
    public class DualSolverLogger
    {
        public int NumCornerDofs { get; set; }
        public int NumExpandedDomainFreeDofs { get; set; }
        public int NumLagrangeMultipliers { get; set; }
        public int NumUniqueGlobalFreeDofs { get; set; }
        public int PcgIterations { get; set; }
        public double PcgResidualNormRatio { get; set; }
    }
}
