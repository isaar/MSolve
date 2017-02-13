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
    class CrackBody2D : IEnrichmentItem2D
    {
        public ICurve2D Geometry { get; }
        public IReadOnlyList<IEnrichmentFunction2D> EnrichmentFunctions { get; }

        public IReadOnlyList<Element2D> AffectedElements { get; private set; }

        public void EnrichNodes()
        {
            // Find all unique affected nodes.
            HashSet<XNode2D> nodes = new HashSet<XNode2D>();
            foreach (var element in AffectedElements) nodes.UnionWith(element.Nodes);
            
            foreach (var node in nodes)
            {
                // Ok for the 1st time. What about updates when only some enrichments must be cleared/changed?
                node.EnrichmentFunctions.UnionWith(this.EnrichmentFunctions); 
            }
        }

    }
}
