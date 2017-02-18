using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Geometry
{
    interface ICurve2D
    {
        double SignedDistanceOf(ICartesianPoint2D point);
        Tuple<double, double> NormalVectorThrough(ICartesianPoint2D point);
    }
}
