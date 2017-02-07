using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry;

namespace ISAAR.MSolve.XFEM.Enrichments
{
    class RampEnrichment2D : IEnrichmentFunction2D
    {
        private readonly ICurve2D discontinuity;

        public RampEnrichment2D(ICurve2D discontinuity)
        {
            this.discontinuity = discontinuity;
        }

        public double ValueAt(IPoint2D point)
        {
            return Math.Abs(discontinuity.SignedDistanceOf(point));
        }

        public Tuple<double, double> DerivativesAt(IPoint2D point)
        {
            int sign = Math.Sign(discontinuity.SignedDistanceOf(point));
            Tuple<double, double> normal = discontinuity.NormalVectorThrough(point);
            return new Tuple<double, double>(sign * normal.Item1, sign * normal.Item2);
        }
    }
}
