using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Functions;

namespace ISAAR.MSolve.XFEM.Entities.FreedomDegrees
{
    class ArtificialDOFType
    {
        public IEnrichmentFunction2D Enrichment { get; }
        public StandardDOFType Axis { get; }

        public ArtificialDOFType(IEnrichmentFunction2D enrichment, StandardDOFType axis)
        {
            this.Enrichment = enrichment;
            this.Axis = axis;
        }
        
        public override string ToString()
        {
            var builder = new StringBuilder();
            builder.Append(Enrichment.ToString());
            builder.Append(' ');
            builder.Append(Axis);
            return builder.ToString();
        }
    }
}
