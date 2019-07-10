using System;
using System.Collections.Generic;
using ISAAR.MSolve.Geometry.Commons;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Geometry.Shapes
{
    public class PolyLine2D : ICurve2D
    {
        private readonly List<double> anglesBetweenSegments;
        private readonly List<double> anglesOfSegments;
        private readonly List<DirectedSegment2D> segments;
        private readonly List<CartesianPoint> vertices;

        public PolyLine2D(CartesianPoint first, CartesianPoint second)
        {
            vertices = new List<CartesianPoint>();
            segments = new List<DirectedSegment2D>();
            anglesBetweenSegments = new List<double>();
            anglesOfSegments = new List<double>();

            vertices.Add(first);
            vertices.Add(second);
            segments.Add(new DirectedSegment2D(first, second));

            double dx = second.X - first.X;
            double dy = second.Y - first.Y;
            anglesOfSegments.Add(Math.Atan2(dy, dx));
        }

        public CartesianPoint End { get { return vertices[vertices.Count - 1]; } }
        public IReadOnlyList<DirectedSegment2D> Segments { get { return segments; } }
        public CartesianPoint Start { get { return vertices[0]; } }
        public Vector2 TangentAtStart => throw new NotImplementedException();
        public Vector2 TangentAtEnd => throw new NotImplementedException();
        public IReadOnlyList<CartesianPoint> Vertices { get { return vertices; } }

        /// <summary>
        /// Counter-clockwise angle from global cartesian x axis to a vector which 1) starts at the end point of the 
        /// curve, 2) is tangent to the curve and 3) heads outwards from the curve.
        /// </summary>
        public double EndPointOrientation() // TODO: perhaps it should return a local coordinate system. I do not need the angle, but its cos and sign.
        {
            CartesianPoint lastSegmentStart = vertices[vertices.Count - 2];
            CartesianPoint lastSegmentEnd = vertices[vertices.Count - 1];
            double dx = lastSegmentEnd.X - lastSegmentStart.X;
            double dy = lastSegmentEnd.Y - lastSegmentStart.Y;
            return Math.Atan2(dy, dx);
        }

        // The normal vector for the positive region.
        public Vector2 NormalVectorThrough(CartesianPoint point)
        {
            if (segments.Count == 1)
            {
                return segments[0].NormalVectorThrough(point);
            }
            else throw new NotImplementedException();
        }

        /// <summary>
        /// See Fries's slides
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        public double SignedDistanceOf(CartesianPoint point)
        {
            if (segments.Count == 1) return segments[0].TransformGlobalToLocalPoint(point).Y;

            var distances = new List<double>();
            bool afterPreviousSegment = false;

            // First segment
            CartesianPoint localPoint = segments[0].TransformGlobalToLocalPoint(point);
            if (localPoint.X < segments[0].Length) distances.Add(localPoint.Y);
            else afterPreviousSegment = true;

            // Subsequent segments
            for (int i = 1; i < segments.Count - 1; ++i)
            {
                localPoint = segments[i].TransformGlobalToLocalPoint(point);
                if (localPoint.X < 0.0)
                {
                    if (afterPreviousSegment)
                    {
                        // Compute the distance from the vertex between this segment and the previous
                        double dx = point.X - vertices[i].X;
                        double dy = point.Y - vertices[i].Y;
                        double distance = Math.Sqrt(dx * dx + dy * dy);
                        int sign = -Math.Sign(anglesBetweenSegments[i - 1]); // If growth angle > 0, the convex angle faces the positive area.
                        distances.Add(sign * distance);
                    }
                    afterPreviousSegment = false;
                }
                else if (localPoint.X <= segments[i].Length)
                {
                    distances.Add(localPoint.Y);
                    afterPreviousSegment = false;
                }
                else afterPreviousSegment = true;
            }

            // Last segment
            int last = segments.Count - 1;
            localPoint = segments[last].TransformGlobalToLocalPoint(point);
            if (localPoint.X < 0.0)
            {
                if (afterPreviousSegment)
                {
                    // Compute the distance from the vertex between this segment and the previous
                    double dx = point.X - vertices[last].X;
                    double dy = point.Y - vertices[last].Y;
                    double distance = Math.Sqrt(dx * dx + dy * dy);
                    int sign = -Math.Sign(anglesBetweenSegments[last - 1]); // If growth angle > 0, the convex angle faces the positive area.
                    distances.Add(sign * distance);
                }
                afterPreviousSegment = false;
            }
            else distances.Add(localPoint.Y);

            return distances[MathUtilities.IndexOfMinAbs(distances)];
        }

        /// <summary>
        /// Counter-clockwise angle from global cartesian x axis to a vector which 1) starts at the start point of the 
        /// curve, 2) is tangent to the curve and 3) heads outwards from the curve.
        /// </summary>
        public double StartPointOrientation()
        {
            // This one's orientation requires more thought, especially since the convention for determining the
            // level set value gives the opposite sign from the rest of the curve.
            CartesianPoint firstSegmentStart = vertices[0];
            CartesianPoint firstSegmentEnd = vertices[1];
            double dx = firstSegmentStart.X - firstSegmentEnd.X;
            double dy = firstSegmentStart.Y - firstSegmentEnd.Y;
            return Math.Atan2(dy, dx);
        }

        public void UpdateGeometry(double angleToLastSegment, double length)
        {
            double lastGlobalAngle = anglesOfSegments[anglesOfSegments.Count - 1];
            double newGlobalAngle = MathUtilities.WrapAngle(angleToLastSegment + lastGlobalAngle);
            double dx = length * Math.Cos(newGlobalAngle);
            double dy = length * Math.Sin(newGlobalAngle);

            var lastPoint = Vertices[Vertices.Count - 1];
            var newPoint = new CartesianPoint(lastPoint.X + dx, lastPoint.Y + dy);
            vertices.Add(newPoint);
            segments.Add(new DirectedSegment2D(lastPoint, newPoint));
            anglesBetweenSegments.Add(angleToLastSegment); // These are independent of the global coordinate system
            anglesOfSegments.Add(newGlobalAngle);
        }

        // Perhaps geometry classes should be decoupled from elements and interact through polygons instead.
        //IReadOnlyList<CartesianPoint> IntersectionWith(XContinuumElement2D element);
        
        // This is the correct one but it needs constrained Delauny triangulation
        //public void IntersectionWith(ContinuumElement2D element, out CartesianPoint[] points, out LineSegment2D[] segments);
    }
}
