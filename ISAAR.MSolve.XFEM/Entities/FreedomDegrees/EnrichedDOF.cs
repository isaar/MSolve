using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Functions;

namespace ISAAR.MSolve.XFEM.Entities.FreedomDegrees
{
    internal class EnrichedDOF: IDOF
    {
        public IEnrichmentFunction2D Enrichment { get; }
        public IDOF StandardDOF { get; }

        public EnrichedDOF(IEnrichmentFunction2D enrichment, IDOF standardDOF)
        {
            this.Enrichment = enrichment;
            this.StandardDOF = standardDOF;
        }
        
        public override string ToString()
        {
            var builder = new StringBuilder();
            builder.Append(Enrichment.ToString());
            builder.Append(" enriched ");
            builder.Append(StandardDOF);
            return builder.ToString();
        }
    }
}
