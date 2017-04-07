using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.Shapes;


namespace ISAAR.MSolve.XFEM.Enrichments.Items.CrackTip
{
    class FixedEnrichmentArea: ITipEnrichmentAreaStrategy
    {
        private readonly double fixedRadius;

        public FixedEnrichmentArea(double enrichmentRadius)
        {
            if (enrichmentRadius <= 0.0) throw new ArgumentException("The radius of the enrichment area must be " 
                + "positive, but was: " + enrichmentRadius);
            this.fixedRadius = enrichmentRadius;
        }

        public double ComputeRadiusOfJintegralOuterContour(CrackTip2D tipItem)
        {
            return fixedRadius;
        }

        public IReadOnlyList<XNode2D> SelectNodesForEnrichment(CrackTip2D tipItem)
        {
            var circle = new Circle2D(tipItem.TipCoordinates, fixedRadius);
            return tipItem.Mesh.FindNodesInsideCircle(circle, true, tipItem.TipElement);
        }
    }
}
