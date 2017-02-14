using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Functions;

namespace ISAAR.MSolve.XFEM.Entities
{
    class XNode2D: Node2D
    {
        public Tuple<IEnrichmentFunction2D, double>[] EnrichmentFunctions { get; set; }
        public bool IsEnriched { get { return EnrichmentFunctions.Length > 0; } }

        // 2 dofs (x,y) for the displacement field and 2 for each enrichment function. 
        // What about for structural elements (rotation) and coupled problems (porous, temperature)?
        public int TotalDofsCount { get { return 2 + ArtificialDofsCount; } }
        public int ArtificialDofsCount { get { return 2 * EnrichmentFunctions.Length; } }


        public XNode2D(int id, double x, double y) : base(id, x, y)
        {
            this.EnrichmentFunctions = new Tuple<IEnrichmentFunction2D, double>[0];
        }
    }
}