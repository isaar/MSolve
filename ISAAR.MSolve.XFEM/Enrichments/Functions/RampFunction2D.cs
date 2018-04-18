using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Geometry;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Enrichments.Functions
{
    class RampFunction2D: IEnrichmentFunction2D
    {
        public RampFunction2D()
        {
        }

        public double EvaluateAt(double signedDistance)
        {
            return Math.Abs(signedDistance);
        }

        public EvaluatedFunction2D EvaluateAllAt(double signedDistance, Tuple<double, double> normalVector)
        {
            int sign = Math.Sign(signedDistance);
            return new EvaluatedFunction2D(Math.Abs(signedDistance), 
                new double[] { sign * normalVector.Item1, sign * normalVector.Item2 });
        }

        public override string ToString()
        {
            return "Ramp";
        }
    }
}
