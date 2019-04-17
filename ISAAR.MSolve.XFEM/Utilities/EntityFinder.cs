using System.Collections.Generic;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.Utilities
{
    class EntityFinder
    {
        private readonly double tolerance; // TODO: Use a value comparer instead

        public EntityFinder(double tolerance = 1e-4)
        {
            this.tolerance = tolerance;
        }
        public List<XContinuumElement2D> FindElementsThatContains(IEnumerable<XContinuumElement2D> elements, 
            CartesianPoint point)
        {
            var result = new List<XContinuumElement2D>();
            foreach (var element in elements)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Nodes);
                PolygonPointPosition position = outline.FindRelativePositionOfPoint(point);
                if (position == PolygonPointPosition.Inside)
                {
                    result.Add(element);
                    break;
                }
                else if ((position == PolygonPointPosition.OnEdge) || (position == PolygonPointPosition.OnVertex))
                {
                    result.Add(element);
                }
            }

            if (result.Count == 0) throw new KeyNotFoundException("No element containing the point " + point + "was found");
            return result;
        }
    }
}
