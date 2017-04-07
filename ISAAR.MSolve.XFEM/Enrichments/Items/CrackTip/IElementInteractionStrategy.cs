using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Enrichments.Items.CrackTip
{
    interface ITipEnrichmentAreaStrategy
    {
        double ComputeRadiusOfJintegralOuterContour(CrackTip2D tipItem);
        IReadOnlyList<XNode2D> SelectNodesForEnrichment(CrackTip2D tipItem);
    }
}
