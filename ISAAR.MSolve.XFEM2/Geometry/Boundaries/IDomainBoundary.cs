using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Geometry.Boundaries
{
    interface IDomainBoundary
    {
        // Not on the boundary exactly.
        bool IsInside(ICartesianPoint2D point);
    }
}
