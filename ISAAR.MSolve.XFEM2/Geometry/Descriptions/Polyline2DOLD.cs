using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.Geometry.Descriptions
{
    class Polyline2DOLD: IGeometryDescription2D
    {
        private readonly List<ICartesianPoint2D> vertices;
        private readonly List<DirectedSegment2D> segments;

        public Dictionary<XContinuumElement2D, CartesianPoint2D[]> ElementIntersections { get; }
        public ICartesianPoint2D StartPoint { get { return vertices[0]; } }
        public ICartesianPoint2D EndPoint { get { return vertices[vertices.Count - 1]; } }

        /// <summary>
        /// Counter-clockwise angle from global cartesian x axis to a vector which 1) starts at the end point of the 
        /// curve, 2) is tangent to the curve and 3) heads outwards from the curve.
        /// </summary>
        public double EndPointOrientation()
        {
            ICartesianPoint2D lastSegmentStart = vertices[vertices.Count - 2];
            ICartesianPoint2D lastSegmentEnd = vertices[vertices.Count - 1];
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
            ICartesianPoint2D firstSegmentStart = vertices[0];
            ICartesianPoint2D firstSegmentEnd = vertices[1];
            double dx = firstSegmentStart.X - firstSegmentEnd.X;
            double dy = firstSegmentStart.Y - firstSegmentEnd.Y;
            return Math.Atan2(dy, dx);
        }

        public Polyline2DOLD(ICartesianPoint2D start, ICartesianPoint2D end)
        {
            vertices = new List<ICartesianPoint2D>();
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
        public double SignedDistanceOf(ICartesianPoint2D point)
        {
            if (segments.Count == 1)
            {
                return segments[0].SignedDistanceOf(point);
            }
            else throw new NotImplementedException();
        }

        // The normal vector for the positive region.
        public Vector2 NormalVectorThrough(ICartesianPoint2D point)
        {
            if (segments.Count == 1)
            {
                return segments[0].NormalVectorThrough(point);
            }
            else throw new NotImplementedException();
        }

        public IReadOnlyList<ICartesianPoint2D> IntersectionWith(XContinuumElement2D element)
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
