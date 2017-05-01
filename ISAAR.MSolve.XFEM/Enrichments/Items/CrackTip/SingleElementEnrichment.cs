using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.Enrichments.Items.CrackTip
{
    class SingleElementEnrichment: ITipEnrichmentAreaStrategy
    {
        public SingleElementEnrichment()
        {
        }

        public IReadOnlyList<XNode2D> SelectNodesForEnrichment(CrackTip2D tipItem)
        {
            var nodes = new List<XNode2D>();
            foreach (var element in tipItem.TipElements) nodes.AddRange(element.Nodes);
            return nodes;
        }
    }
}
