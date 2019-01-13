using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Enrichments.Functions
{
    interface IHeavisideFunction2D : IEnrichmentFunction2D
    {
        double EvaluateAt(double signedDistance);
        EvaluatedFunction2D EvaluateAllAt(double signedDistance);
    }
}
