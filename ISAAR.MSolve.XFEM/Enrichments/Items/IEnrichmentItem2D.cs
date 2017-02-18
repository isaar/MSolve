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
    // Connects the geometry, model and enrichment function entities.
    interface IEnrichmentItem2D
    {
        IReadOnlyList<IEnrichmentFunction2D> EnrichmentFunctions { get; }
        IReadOnlyList<XElement2D> AffectedElements { get; }

        void AffectElement(XElement2D element);

        /// <summary>
        /// Assigns enrichment functions and their nodal values to each enriched node.
        /// </summary>
        void EnrichNodes();
    }
}
