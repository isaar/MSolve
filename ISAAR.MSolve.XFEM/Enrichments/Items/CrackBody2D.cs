using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Geometry;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    class CrackBody2D : AbstractEnrichmentItem2D
    {
        public ICurve2D Discontinuity { get; }

        //TODO: allow the user to specify which Heaviside function to pass. 
        // First specialize IEnrichmentFunction2D into IDiscontinuityFunction2D to allow only Heaviside, Sign or 
        // Heaviside approximation functions.
        public CrackBody2D(ICurve2D discontinuity)
        {
            this.Discontinuity = discontinuity;
            IEnrichmentFunction2D enrichment = new SignFunction2D(this);
            this.EnrichmentFunctions = new IEnrichmentFunction2D[] { enrichment };
            this.DOFs = new ArtificialDOFType[] {
                new ArtificialDOFType(enrichment, StandardDOFType.X),
                new ArtificialDOFType(enrichment, StandardDOFType.Y)
            };
        }

        public override IReadOnlyList<ICartesianPoint2D> IntersectionPointsForIntegration(XElement2D element)
        {
            return Discontinuity.IntersectionWith(element);
        }
    }
}
