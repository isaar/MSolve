using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry;

namespace ISAAR.MSolve.XFEM.Enrichments
{
    class HeavisideEnrichment2D : IEnrichmentFunction2D
    {
        private readonly ICurve2D discontinuity;

        public HeavisideEnrichment2D(ICurve2D discontinuity)
        {
            this.discontinuity = discontinuity;
        }

        public double ValueAt(IPoint2D point)
        {
            return (discontinuity.SignedDistanceOf(point) >= 0) ? 1.0 : -1.0; // What should H(0) be? For now H(0) = 1.
        }

        public Tuple<double, double> DerivativesAt(IPoint2D point)
        {
            return new Tuple<double, double>(0.0, 0.0);
        }
    }
}
