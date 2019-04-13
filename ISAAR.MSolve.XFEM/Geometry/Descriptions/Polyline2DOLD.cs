using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.Geometry.Descriptions
{
    class Polyline2DOLD: IGeometryDescription2D
    {
        private readonly List<CartesianPoint2D> vertices;
        private readonly List<DirectedSegment2D> segments;

        public Dictionary<XContinuumElement2D, CartesianPoint2D[]> ElementIntersections { get; }
        public CartesianPoint2D StartPoint { get { return vertices[0]; } }
        public CartesianPoint2D EndPoint { get { return vertices[vertices.Count - 1]; } }

        /// <summary>
        /// Counter-clockwise angle from global cartesian x axis to a vector which 1) starts at the end point of the 
        /// curve, 2) is tangent to the curve and 3) heads outwards from the curve.
        /// </summary>
        public double EndPointOrientation()
        {
            CartesianPoint2D lastSegmentStart = vertices[vertices.Count - 2];
            CartesianPoint2D lastSegmentEnd = vertices[vertices.Count - 1];
            double dx = lastSegmentEnd.X - lastSegmentStart.X;
            double dy = lastSegmentEnd.Y - lastSegmentStart.Y;
            return Math.Atan2(dy, dx);
        }

        /// <summary>
        /// Counter-clockwise angle from global cartesian x axis to a vector which 1) starts at the start point of the 
        /// curve, 2) is tangent to the curve and 3) heads outwards from the curve.
        /// </summary>
        public double StartPointOrientation()
        {
            CartesianPoint2D firstSegmentStart = vertices[0];
            CartesianPoint2D firstSegmentEnd = vertices[1];
            double dx = firstSegmentStart.X - firstSegmentEnd.X;
            double dy = firstSegmentStart.Y - firstSegmentEnd.Y;
            return Math.Atan2(dy, dx);
        }

        public Polyline2DOLD(CartesianPoint2D start, CartesianPoint2D end)
        {
            vertices = new List<CartesianPoint2D>();
            vertices.Add(start);
            vertices.Add(end);

            segments = new List<DirectedSegment2D>();
            segments.Add(new DirectedSegment2D(start, end));

            ElementIntersections = new Dictionary<XContinuumElement2D, CartesianPoint2D[]>();
        }

        /// <summary>
        /// See Fries's slides
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        public double SignedDistanceOf(CartesianPoint2D point)
        {
            if (segments.Count == 1)
            {
                return segments[0].SignedDistanceOf(point);
            }
            else throw new NotImplementedException();
        }

        // The normal vector for the positive region.
        public Vector2 NormalVectorThrough(CartesianPoint2D point)
        {
            if (segments.Count == 1)
            {
                return segments[0].NormalVectorThrough(point);
            }
            else throw new NotImplementedException();
        }

        public IReadOnlyList<CartesianPoint2D> IntersectionWith(XContinuumElement2D element)
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
