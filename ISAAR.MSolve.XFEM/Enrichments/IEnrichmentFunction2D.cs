using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry;

namespace ISAAR.MSolve.XFEM.Enrichments
{
    interface IEnrichmentFunction2D
    {
        double ValueAt(IPoint2D point);
        Tuple<double, double> DerivativesAt(IPoint2D point);
    }
}
