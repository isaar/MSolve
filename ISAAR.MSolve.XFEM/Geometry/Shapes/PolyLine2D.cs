using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Utilities;

//TODO: implement all inherited methods
//TODO: remove dependencies from classes developed for crack propagation
namespace ISAAR.MSolve.XFEM.Geometry.Shapes
{
    class PolyLine2D : IOpenCurve2D
    {
        private readonly List<double> angles;
        private readonly List<DirectedSegment2D> segments;
        private readonly List<ICartesianPoint2D> vertices;
        private TipCoordinateSystem tipSystem;

        public PolyLine2D(ICartesianPoint2D first, ICartesianPoint2D second)
        {
            vertices = new List<ICartesianPoint2D>();
            segments = new List<DirectedSegment2D>();
            angles = new List<double>();

            double dx = second.X - first.X;
            double dy = second.Y - first.Y;
            double tangentSlope = Math.Atan2(dy, dx);
            tipSystem = new TipCoordinateSystem(second, tangentSlope);

            vertices.Add(first);
            vertices.Add(second);
            segments.Add(new DirectedSegment2D(first, second));
        }

        public IReadOnlyList<DirectedSegment2D> Segments { get { return segments; } }
        public IReadOnlyList<ICartesianPoint2D> Vertices { get { return vertices; } }
        public ICartesianPoint2D Start { get { return vertices[0]; } }
        public ICartesianPoint2D End { get { return vertices[vertices.Count - 1]; } }
        public Vector2 TangentAtStart => throw new NotImplementedException();
        public Vector2 TangentAtEnd => throw new NotImplementedException();

        // The normal vector for the positive region.
        public Vector2 NormalThrough(ICartesianPoint2D point)
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
        public double SignedDistance(ICartesianPoint2D point)
        {
            if (segments.Count == 1) return segments[0].TransformGlobalToLocalPoint(point).Y;

            var distances = new List<double>();
            bool afterPreviousSegment = false;

            // First segment
            ICartesianPoint2D localPoint = segments[0].TransformGlobalToLocalPoint(point);
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
                        int sign = -Math.Sign(angles[i - 1]); // If growth angle > 0, the convex angle faces the positive area.
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
                    int sign = -Math.Sign(angles[last - 1]); // If growth angle > 0, the convex angle faces the positive area.
                    distances.Add(sign * distance);
                }
                afterPreviousSegment = false;
            }
            else distances.Add(localPoint.Y);

            return distances[MathUtilities.IndexOfMinAbs(distances)];
        }

        public void UpdateGeometry(double localGrowthAngle, double growthLength)
        {
            double globalGrowthAngle = AngleUtilities.Wrap(localGrowthAngle + tipSystem.RotationAngle);
            double dx = growthLength * Math.Cos(globalGrowthAngle);
            double dy = growthLength * Math.Sin(globalGrowthAngle);

            var oldTip = Vertices[Vertices.Count - 1];
            var newTip = new CartesianPoint2D(oldTip.X + dx, oldTip.Y + dy);
            vertices.Add(newTip);
            segments.Add(new DirectedSegment2D(oldTip, newTip));
            angles.Add(localGrowthAngle); // These are independent of the global coordinate system
            tipSystem = new TipCoordinateSystem(newTip, globalGrowthAngle);
        }
    }
}
