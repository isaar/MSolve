using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Geometry;
using ISAAR.MSolve.XFEM.Utilities;

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
        private readonly MaterialInterface2D enrichmentItem;

        public RampFunction2D(MaterialInterface2D enrichmentItem)
        {
            this.enrichmentItem = enrichmentItem;
        }

        public double EvalueAt(IPoint2D cartesianPoint)
        {
            return Math.Abs(enrichmentItem.Discontinuity.SignedDistanceOf(cartesianPoint));
        }

        public Tuple<double, double> EvaluateDerivativesAt(IPoint2D cartesianPoint)
        {
            int sign = Math.Sign(enrichmentItem.Discontinuity.SignedDistanceOf(cartesianPoint));
            Tuple<double, double> normal = enrichmentItem.Discontinuity.NormalVectorThrough(cartesianPoint);
            return new Tuple<double, double>(sign * normal.Item1, sign * normal.Item2);
        }

        public EvaluatedFunction2D EvaluateAllAt(IPoint2D cartesianPoint)
        {
            double signedDistance = enrichmentItem.Discontinuity.SignedDistanceOf(cartesianPoint);
            int sign = Math.Sign(signedDistance);
            Tuple<double, double> normal = enrichmentItem.Discontinuity.NormalVectorThrough(cartesianPoint);
            return new EvaluatedFunction2D(Math.Abs(signedDistance), 
                new Tuple<double, double>(sign * normal.Item1, sign * normal.Item2));
        }
    }
}
