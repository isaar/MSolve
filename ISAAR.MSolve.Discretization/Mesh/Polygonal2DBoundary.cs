using System.Collections.Generic;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;

namespace ISAAR.MSolve.Discretization.Mesh
{
    public class Polygonal2DBoundary: IDomain2DBoundary
    {
        private readonly ConvexPolygon2D polygon;

        public Polygonal2DBoundary(IReadOnlyList<CartesianPoint> vertices)
        {
            polygon = ConvexPolygon2D.CreateUnsafe(vertices); 
        }

        public bool IsInside(CartesianPoint point)
        {
            var pos = polygon.FindRelativePositionOfPoint(point);
            if (pos == PolygonPointPosition.Inside) return true;
            else return false;
        }
    }
}
