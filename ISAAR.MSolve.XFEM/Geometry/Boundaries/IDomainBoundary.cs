using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.XFEM.Geometry.Boundaries
{
    interface IDomainBoundary
    {
        // Not on the boundary exactly.
        bool IsInside(CartesianPoint2D point);
    }
}
