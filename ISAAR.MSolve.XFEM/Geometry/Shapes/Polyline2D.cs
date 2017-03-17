using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Geometry.Shapes
{
    class Polyline2D: ICurve2D
    {
        private readonly List<ICartesianPoint2D> vertices;
        private readonly List<LineSegment2D> segments;

        public Dictionary<XElement2D, CartesianPoint2D[]> ElementIntersections { get; }

        public Polyline2D(ICartesianPoint2D start, ICartesianPoint2D end)
        {
            vertices = new List<ICartesianPoint2D>();
            vertices.Add(start);
            vertices.Add(end);

            segments = new List<LineSegment2D>();
            segments.Add(new LineSegment2D(start, end));

            ElementIntersections = new Dictionary<XElement2D, CartesianPoint2D[]>();
        }

        /// <summary>
        /// See Fries's slides
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        public double SignedDistanceOf(ICartesianPoint2D point)
        {
           throw new NotImplementedException();
        }

        // The normal vector for the positive region.
        public Tuple<double, double> NormalVectorThrough(ICartesianPoint2D point)
        {
            throw new NotImplementedException();
        }

        public void IntersectionWith(ContinuumElement2D element, out ICartesianPoint2D[] points, out LineSegment2D[] segments)
        {
            throw new NotImplementedException();
        }

        public IReadOnlyList<ICartesianPoint2D> IntersectionWith(XElement2D element)
        {
            CartesianPoint2D[] intersectionPoints;
            bool alreadyIntersected = ElementIntersections.TryGetValue(element, out intersectionPoints);
            if (alreadyIntersected) return intersectionPoints;
            else
            {
                throw new NotImplementedException();
            }
        }
    }
}
