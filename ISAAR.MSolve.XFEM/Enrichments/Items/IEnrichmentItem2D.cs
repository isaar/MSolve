using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    // Connects the geometry, model and enrichment function entities.
    interface IEnrichmentItem2D
    {
        ICurve2D Geometry { get; }
        IReadOnlyList<IEnrichmentFunction2D> EnrichmentFunctions { get; }
        IReadOnlyList<Element2D> AffectedElements { get; }
        void EnrichNodes();
    }
}
