using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    class MaterialInterface2D : AbstractEnrichmentItem2D
    {
        public ICurve2D Discontinuity { get; }

        public MaterialInterface2D(ICurve2D geometry)
        {
            this.Discontinuity = geometry;
            this.EnrichmentFunctions = new IEnrichmentFunction2D[] { new RampFunction2D(this) };
        }
    }
}
