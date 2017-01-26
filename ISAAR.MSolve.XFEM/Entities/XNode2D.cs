using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments;

namespace ISAAR.MSolve.XFEM.Entities
{
    class XNode2D: Node2D
    {
        public ISet<IEnrichmentFunction2D> Enrichments { get; }
        public bool IsEnriched { get { return Enrichments.Count > 0; } }

        public XNode2D(int id, double x, double y): base(id, x, y)
        {
            this.Enrichments = new HashSet<IEnrichmentFunction2D>();
        }
    }
}
