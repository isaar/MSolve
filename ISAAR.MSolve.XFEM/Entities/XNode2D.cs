using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Entities
{
    class XNode2D: Node2D
    {
        public List<Tuple<IEnrichmentFunction2D, double>> EnrichmentFunctions { get; set; }
        public List<IEnrichmentItem2D> EnrichmentItems { get; set; }
        public bool IsEnriched { get { return EnrichmentFunctions.Count > 0; } }

        public int ArtificialDofsCount { get { return 2 * EnrichmentFunctions.Count; } }

        public XNode2D(int id, double x, double y) : base(id, x, y)
        {
            this.EnrichmentFunctions = new List<Tuple<IEnrichmentFunction2D, double>>();
            this.EnrichmentItems = new List<IEnrichmentItem2D>();
        }
    }
}