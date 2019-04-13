using System.Collections.Generic;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.Geometry.Boundaries
{
    class PolygonalBoundary: IDomainBoundary
    {
        private readonly ConvexPolygon2D polygon;

        public PolygonalBoundary(IReadOnlyList<CartesianPoint2D> vertices)
        {
            polygon = ConvexPolygon2D.CreateUnsafe(vertices); 
        }

        public bool IsInside(CartesianPoint2D point)
        {
            var pos = polygon.FindRelativePositionOfPoint(point);
            if (pos == PolygonPointPosition.Inside) return true;
            else return false;
        }
    }
}
