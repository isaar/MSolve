using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Enrichments.Functions
{
    interface ITipFunction: IEnrichmentFunction2D
    {
        double EvaluateAt(PolarPoint2D point);
        EvaluatedFunction2D EvaluateAllAt(PolarPoint2D point, TipJacobians jacobian);
    }
}
