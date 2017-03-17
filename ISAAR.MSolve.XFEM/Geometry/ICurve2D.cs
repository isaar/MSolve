using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.Geometry
{
    interface ICurve2D
    {
        double SignedDistanceOf(ICartesianPoint2D point);
        Tuple<double, double> NormalVectorThrough(ICartesianPoint2D point);

        // Perhaps geometry classes should be decoupled from elements and interact through polygons instead.
        IReadOnlyList<ICartesianPoint2D> IntersectionWith(XElement2D element);

        // This is the correct one but it needs constrained Delauny triangulation
        //public void IntersectionWith(ContinuumElement2D element, out ICartesianPoint2D[] points, out LineSegment2D[] segments);
    }
}
