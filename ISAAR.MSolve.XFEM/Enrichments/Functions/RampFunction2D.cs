using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Geometry;

namespace ISAAR.MSolve.XFEM.Enrichments.Functions
{
    // Cons: 
    // 1) Either I store the enrichment item or the geometry entity. Either way this class (and any similar) is coupled
    //      with geometry.
    // 2) Evaluating both the derivatives and the raw value requires 2 identical signed distance calculations, which 
    //      may be costly. Also the normal vector through a point and that point's signed distance may have common 
    //      costly calculations. It seems that the curve should somehow remember (and clear) the various geometric 
    //      values associated each point.
    class RampFunction2D: IEnrichmentFunction2D
    {
        private readonly IEnrichmentItem2D enrichmentItem;

        public RampFunction2D(IEnrichmentItem2D enrichmentItem)
        {
            this.enrichmentItem = enrichmentItem;
        }

        public double ValueAt(IPoint2D point)
        {
            return Math.Abs(enrichmentItem.Geometry.SignedDistanceOf(point));
        }

        public Tuple<double, double> DerivativesAt(IPoint2D point)
        {
            int sign = Math.Sign(enrichmentItem.Geometry.SignedDistanceOf(point));
            Tuple<double, double> normal = enrichmentItem.Geometry.NormalVectorThrough(point);
            return new Tuple<double, double>(sign * normal.Item1, sign * normal.Item2);
        }
    }
}
