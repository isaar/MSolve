using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    class PropagationResults
    {
        public double crackAngleDegrees;
        public double tipElementSize;
        public double jIntegralRadiusOverElementSize;
        public double enrichmentRadiusOverElementSize;
        public double? startJIntegral;
        public double? startSIF1;
        public double? startSIF2;
        public double endJIntegral;
        public double endSIF1;
        public double endSIF2;
        public double referenceJIntegral;
        public double referenceSIF1;
        public double referenceSIF2;
    }
}
