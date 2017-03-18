﻿using System;
using System.Collections.Generic;
using System.Windows;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Geometry.Shapes
{
    // The line divides the 2D space into a positive and a negative region. Points in the positive region have positive
    // signed distance to the line, while points in the negative region have negative signed distance. Each region also
    // has a normal vector. The normal of the positive region is defined by the following CONVENTION:
    // In a local 1D coordinate system, the START point of the line is to the left, the END point is to the right and
    // the normal vector of the positive region is pointing DOWNWARDS (into the positive region). 
    class LineSegment2D : ICurve2D
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

        public double SignedDistanceOf(ICartesianPoint2D point)
        {
            // The following area is positive in the positive region. TODO: prove it
            double triangleAreax2 = 
                (End.Y - Start.Y) * point.X - (End.X - Start.X) * point.Y + End.X * Start.Y - End.Y * Start.X;
            return triangleAreax2 / Length;
        }

        // The normal vector for the positive region.
        public Tuple<double, double> NormalVectorThrough(ICartesianPoint2D point)
        {
            double dy = End.Y - Start.Y;
            double dx = Start.X - End.X;
            double length = Math.Sqrt(dx * dx + dy * dy);
            return new Tuple<double, double>(dy / length, -dx / length);
        }

        public IReadOnlyList<ICartesianPoint2D> IntersectionWith(XContinuumElement2D element)
        {
            var intersectionPoints = new List<CartesianPoint2D>();
            // Should I also include the vertices if they fall inside? Or should I do that in the enrichment item?
            var polygon = ConvexPolygon2D.CreateUnsafe(element.Nodes);
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
            Vector p = new Vector(this.Start.X, this.Start.Y);
            Vector r = new Vector(this.End.X - this.Start.X, this.End.Y - this.Start.Y);
            Vector q = new Vector(segment.Start.X, segment.Start.Y);
            Vector s = new Vector(segment.End.X - segment.Start.X, segment.End.Y - segment.Start.Y);

            // Find the cross products r x s and (q-p)*r 
            double rCrossS = Vector.CrossProduct(r, s);
            Vector qMinusP = Vector.Subtract(q, p);
            double qMinusPCrossR = Vector.CrossProduct(qMinusP, r);

            if (rCrossS == 0) // TODO: use tolerance here
            {
                if (qMinusPCrossR == 0) // TODO: use tolerance here
                {
                    double rDotR = r.LengthSquared;
                    double t0 = Vector.Multiply(qMinusP, r) / rDotR;
                    double t1 = t0 + Vector.Multiply(s, r) / rDotR;
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
                double t = Vector.CrossProduct(qMinusP, s) / rCrossS;
                double u = qMinusPCrossR / rCrossS;
                if ((t >= 0.0) && (t <= 1.0) && (u >= 0.0) && (u <= 1.0))
                {
                    Vector solution = Vector.Add(p, Vector.Multiply(t, r));
                    intersectionPoint = new CartesianPoint2D(solution.X, solution.Y);
                    return SegmentSegmentPosition.Intersecting;
                }
                else
                {
                    intersectionPoint = null;
                    return SegmentSegmentPosition.Disjoint;
                }
            }
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
