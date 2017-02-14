using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Enrichments.Functions
{
    // TODO: Purge the derivatives only evaluation method, when no longer needed.
    interface IEnrichmentFunction2D
    {
        double EvalueAt(IPoint2D cartesianPoint);
        Tuple<double, double> EvaluateDerivativesAt(IPoint2D cartesianPoint);
        EvaluatedFunction2D EvaluateAllAt(IPoint2D cartesianPoint);
    }
}
