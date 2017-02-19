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
    // Cons: 
    // 1) Either I store the enrichment item or the geometry entity. Either way this class (and any similar) is coupled
    //      with geometry.
    // 2) Evaluating both the derivatives and the raw value requires 2 identical signed distance calculations, which 
    //      may be costly. Also the normal vector through a point and that point's signed distance may have common 
    //      costly calculations. It seems that the curve should somehow remember (and clear) the various geometric 
    //      values associated each point.
    class SignFunction2D : IEnrichmentFunction2D
    {
        private readonly CrackBody2D enrichmentItem;

        public SignFunction2D(CrackBody2D enrichmentItem)
        {
            this.enrichmentItem = enrichmentItem;
        }

        public double EvalueAt(ICartesianPoint2D cartesianPoint)
        {
            double signedDistance = enrichmentItem.Discontinuity.SignedDistanceOf(cartesianPoint);
            if (signedDistance > 0.0) return 1.0;
            else if (signedDistance < 0.0) return -1.0;
            else return 0.0;
        }

        public EvaluatedFunction2D EvaluateAllAt(ICartesianPoint2D cartesianPoint)
        {
            var derivatives = new Tuple<double, double>(0.0, 0.0);
            double signedDistance = enrichmentItem.Discontinuity.SignedDistanceOf(cartesianPoint);
            if (signedDistance > 0) return new EvaluatedFunction2D(1.0, derivatives);
            else if (signedDistance < 0) return new EvaluatedFunction2D(-1.0, derivatives);
            else return new EvaluatedFunction2D(0.0, derivatives);
        }
    }
}
