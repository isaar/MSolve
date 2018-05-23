using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.Geometry.Boundaries
{
    class PolygonalBoundary: IDomainBoundary
    {
        private readonly ConvexPolygon2D polygon;

        public PolygonalBoundary(IReadOnlyList<ICartesianPoint2D> vertices)
        {
            polygon = ConvexPolygon2D.CreateUnsafe(vertices); 
        }

        public bool IsInside(ICartesianPoint2D point)
        {
            var pos = polygon.FindRelativePositionOfPoint(point);
            if (pos == PolygonPointPosition.Inside) return true;
            else return false;
        }
    }
}
