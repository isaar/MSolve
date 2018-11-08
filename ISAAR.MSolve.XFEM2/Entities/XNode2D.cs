using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Entities
{
    class XNode2D: Node2D
    {
        public Dictionary<IEnrichmentItem2D, double[]> EnrichmentItems { get; }
        public bool IsEnriched { get { return EnrichmentItems.Count > 0; } }

        public int EnrichedDofsCount
        {
            get
            {
                int count = 0;
                foreach (IEnrichmentItem2D enrichment in EnrichmentItems.Keys)
                {
                    count += enrichment.Dofs.Count;
                }
                return count;
            }
        }

        public ISet<EnrichedDof> EnrichedDofs
        {
            get
            {
                var dofs = new HashSet<EnrichedDof>();
                foreach (IEnrichmentItem2D enrichment in EnrichmentItems.Keys)
                {
                    dofs.UnionWith(enrichment.Dofs);
                }
                return dofs;
            }
        }

        public XNode2D(int id, double x, double y) : base(id, x, y)
        {
            this.EnrichmentItems = new Dictionary<IEnrichmentItem2D, double[]>();
        }
    }
}