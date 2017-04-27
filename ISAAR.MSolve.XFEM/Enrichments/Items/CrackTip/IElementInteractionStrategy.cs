using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Enrichments.Items.CrackTip
{
    // TODO: Perhaps I should decouple the J-integral radius form the enrichment radius. There is no real theoritical 
    // reason to have them be the same.
    interface ITipEnrichmentAreaStrategy
    {
        double ComputeRadiusOfJintegralOuterContour(CrackTip2D tipItem);
        IReadOnlyList<XNode2D> SelectNodesForEnrichment(CrackTip2D tipItem);
    }
}
