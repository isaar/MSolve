using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.XFEM.Enrichments.Functions;

namespace ISAAR.MSolve.XFEM.FreedomDegrees
{
    internal class EnrichedDof: IDof
    {
        public IEnrichmentFunction2D Enrichment { get; }
        public IDof StandardDof { get; }

        public EnrichedDof(IEnrichmentFunction2D enrichment, IDof standardDof)
        {
            this.Enrichment = enrichment;
            this.StandardDof = standardDof;
        }
        
        public override string ToString()
        {
            var builder = new StringBuilder();
            builder.Append(Enrichment.ToString());
            builder.Append(" enriched ");
            builder.Append(StandardDof);
            return builder.ToString();
        }
    }
}
