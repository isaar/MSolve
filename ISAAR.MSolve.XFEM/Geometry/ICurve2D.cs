using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Geometry
{
    interface ICurve2D
    {
        double SignedDistanceOf(IPoint2D point);
    }
}
