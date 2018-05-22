using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Geometry.Shapes
{
    // The line divides the 2D space into a positive and a negative region. Points in the positive region have positive
    // signed distance to the line, while points in the negative region have negative signed distance. Each region also
    // has a normal vector. The normal of the positive region is defined by the following CONVENTION:
    // In a local 1D coordinate system, the START point of the line is to the left, the END point is to the right and
    // the normal vector of the positive region is pointing DOWNWARDS (into the positive region). 
    class LineSegment2D
    {
        public ICartesianPoint2D Start { get; }
        public ICartesianPoint2D End { get; }

        public double Length
        {
            get
            {
                double dx = End.X - Start.X;
                double dy = End.Y - Start.Y;
                return Math.Sqrt(dx * dx + dy * dy);
            }
        }
        
        public LineSegment2D(ICartesianPoint2D start, ICartesianPoint2D end)
        {
            this.Start = start;
            this.End = end;
        }

        public double DistanceOf(ICartesianPoint2D point)
        {
            double triangleAreax2 = Math.Abs(
                (End.Y - Start.Y) * point.X - (End.X - Start.X) * point.Y + End.X * Start.Y - End.Y * Start.X);
            return triangleAreax2 / Length;
        }

        // One of the 2 normal vectors for the positive region.
        public Vector NormalVectorThrough(ICartesianPoint2D point)
        {
            double dy = End.Y - Start.Y;
            double dx = Start.X - End.X;
            double length = Math.Sqrt(dx * dx + dy * dy);
            return Vector.CreateFromArray(new double[] { dy / length, -dx / length });
        }

        public IReadOnlyList<ICartesianPoint2D> IntersectionWith(ConvexPolygon2D polygon)
        {
            var intersectionPoints = new List<CartesianPoint2D>();
            // Should I also include the vertices if they fall inside? Or should I do that in the enrichment item?
            foreach (LineSegment2D edge in polygon.Edges)
            {
                CartesianPoint2D point;
                SegmentSegmentPosition intersection = IntersectionWith(edge, out point);
                if (intersection == SegmentSegmentPosition.Intersecting) intersectionPoints.Add(point);
            }
            return intersectionPoints;
        }

        /// <summary>
        /// See <see cref="http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect"/>
        /// </summary>
        /// <param name="segment"></param>
        /// <param name="intersectionPoint"></param>
        /// <returns></returns>
        public SegmentSegmentPosition IntersectionWith(LineSegment2D segment, out CartesianPoint2D intersectionPoint) 
        {
            //TODO: optimize this:
            // a) Cache the vector representation of "this".Perhaps there could be a private vector representation class
            //b) Perhaps some cases can be ignored/lumped together.E.g. for colinear
            //c) Perhaps working with the double values instead of vector structs saves a lot of accesses.

            //TODO: handle extreme cases:
            // a) One segment's end lies on the other segment
            // b) One segment's end coincides with one of the other segment's ends, but otherwise they do not overlap
            // c) These cases probably need some tolerance checks.This requires attention since intersection points
            //      outside the element will create unwanted triangulations.

            //Create vectors p, r, q, s
            var p = Vector2.Create(this.Start.X, this.Start.Y);
            var r = Vector2.Create(this.End.X - this.Start.X, this.End.Y - this.Start.Y);
            var q = Vector2.Create(segment.Start.X, segment.Start.Y);
            var s = Vector2.Create(segment.End.X - segment.Start.X, segment.End.Y - segment.Start.Y);

            // Find the cross products r x s and (q-p)*r 
            double rCrossS = r.CrossProduct(s);
            var qMinusP = q - p;
            double qMinusPCrossR = qMinusP.CrossProduct(r);

            if (rCrossS == 0) // TODO: use tolerance here
            {
                if (qMinusPCrossR == 0) // TODO: use tolerance here
                {
                    double rDotR = r * r;
                    double t0 = (qMinusP * r) / rDotR;
                    double t1 = t0 + (s * r) / rDotR;
                    LineSegment1D this1D = (new LineSegment1D(0.0, 1.0));
                    LineSegment1D.SegmentSegmentPosition intersection1D = 
                        this1D.IntesectionWith(new LineSegment1D(t0, t1));
                    if (intersection1D == LineSegment1D.SegmentSegmentPosition.Disjoint) 
                    {
                        intersectionPoint = null;
                        return LineSegment2D.SegmentSegmentPosition.CollinearDisjoint;
                    }
                    else if (intersection1D == LineSegment1D.SegmentSegmentPosition.Overlapping)
                    {
                        intersectionPoint = null;
                        return LineSegment2D.SegmentSegmentPosition.Overlapping;
                    }
                    else // TODO: handle this case
                    {
                        throw new NotImplementedException("The two segments are colinear and meet only at one point");
                    }
                }
                else
                {
                    intersectionPoint = null;
                    return SegmentSegmentPosition.Parallel;
                }
            }
            else
            {
                double t = qMinusP.CrossProduct(s) / rCrossS;
                double u = qMinusPCrossR / rCrossS;
                if ((t >= 0.0) && (t <= 1.0) && (u >= 0.0) && (u <= 1.0))
                {
                    var solution = p + t * r;
                    intersectionPoint = new CartesianPoint2D(solution);
                    return SegmentSegmentPosition.Intersecting;
                }
                else
                {
                    intersectionPoint = null;
                    return SegmentSegmentPosition.Disjoint;
                }
            }
        }

        public SegmentPointPosition FindRelativePositionOfPoint(ICartesianPoint2D point)
        {
            // Local coordinate system
            double length = Length;
            double cosa = (End.X - Start.X) / length;
            double sina = (End.Y - Start.Y) / length;
            double originLocalX = -cosa * Start.X - sina * Start.Y;
            double originLocalY = sina * Start.X - cosa * Start.Y;

            double localX = cosa * point.X + sina * point.Y + originLocalX;
            double localY = -sina * point.X + cosa * point.Y + originLocalY;
            double tolerance = 1e-6;

            if (Math.Abs(localY) < tolerance)
            {
                if ((localX >= -tolerance) && (localX <= length + tolerance)) return SegmentPointPosition.PointOnSegment;
            }
            return SegmentPointPosition.Disjoint;
        }

        public enum SegmentPointPosition
        {
            PointOnSegment, Disjoint
        }

        public enum SegmentSegmentPosition // TODO: Perhaps this should be a class
        {
            Disjoint,
            CollinearDisjoint,
            Parallel,

            /// <summary>
            /// This includes the cases where one segments ends on the other.
            /// </summary>
            Intersecting,

            Overlapping
        }
    }
}
