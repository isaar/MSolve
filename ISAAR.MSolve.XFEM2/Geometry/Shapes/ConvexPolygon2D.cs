using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Geometry.Shapes
{
    public enum CirclePolygonPosition
    {
        Disjoint, Intersecting, PolygonInsideCircle, CircleInsidePolygon
    }

    public enum PolygonPointPosition
    {
        Inside, Outside, OnEdge, OnVertex
    }

    class ConvexPolygon2D
    {
        private ConvexPolygon2D(IReadOnlyList<ICartesianPoint2D> vertices)
        {
            this.Vertices = vertices;

            var edges = new LineSegment2D[vertices.Count];
            for (int i = 0; i < vertices.Count; ++i)
            {
                edges[i] = new LineSegment2D(vertices[i], vertices[(i + 1) % vertices.Count]);
            }
            this.Edges = edges;
        }

        public IReadOnlyList<ICartesianPoint2D> Vertices { get; }
        public IReadOnlyList<LineSegment2D> Edges { get; }

        public static ConvexPolygon2D CreateUnsafe(IReadOnlyList<ICartesianPoint2D> vertices)
        {
            return new ConvexPolygon2D(vertices);
        }

        public static ConvexPolygon2D CreateSafe(IReadOnlyList<ICartesianPoint2D> vertices)
        {
            var points = new CartesianPoint2D[vertices.Count];
            for (int i = 0; i < vertices.Count; ++i)
            {
                points[i] = new CartesianPoint2D(vertices[i].X, vertices[i].Y);
            }
            return new ConvexPolygon2D(points);
        }

        public double ComputeArea()
        {
            double sum = 0.0;
            for (int vertexIdx = 0; vertexIdx < Vertices.Count; ++vertexIdx)
            {
                ICartesianPoint2D vertex1 = Vertices[vertexIdx];
                ICartesianPoint2D vertex2 = Vertices[(vertexIdx + 1) % Vertices.Count];
                sum += vertex1.X * vertex2.Y - vertex2.X * vertex1.Y;
            }
            return Math.Abs(0.5 * sum); // area would be negative if vertices were in counter-clockwise order
        }

        /// <summary>
        /// Ray casting method, shamelessly copied from 
        /// http://stackoverflow.com/questions/8721406/how-to-determine-if-a-point-is-inside-a-2d-convex-polygon or 
        /// http://stackoverflow.com/questions/217578/how-can-i-determine-whether-a-2d-point-is-within-a-polygon?noredirect=1&lq=1
        /// For points on the outline of the polygon, false is returned. NOPE, it isn't always.
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        public bool IsPointInsidePolygon(ICartesianPoint2D point)
        {
            int i, j;
            bool result = false;
            for (i = 0, j = Vertices.Count - 1; i < Vertices.Count; j = i++)
            {
                if (    ( (Vertices[i].Y > point.Y) != (Vertices[j].Y > point.Y) ) &&
                        ( point.X < (Vertices[j].X - Vertices[i].X) * (point.Y - Vertices[i].Y) 
                                    / (Vertices[j].Y - Vertices[i].Y) + Vertices[i].X ) )
                {
                    result = !result;
                }
            }
            return result;
        }

        public CirclePolygonPosition FindRelativePositionOfCircle(Circle2D circle)
        {
            int verticesOutsideCircle = 0;
            int verticesInsideCircle = 0;
            foreach (ICartesianPoint2D vertex in Vertices)
            {
                CirclePointPosition vertexPosition = circle.FindRelativePositionOfPoint(vertex);
                if (vertexPosition == CirclePointPosition.Outside) ++verticesOutsideCircle;
                else if (vertexPosition == CirclePointPosition.Inside) ++verticesInsideCircle;
            }
            int verticesOnCircle = Vertices.Count - verticesOutsideCircle - verticesInsideCircle;

            if (verticesOutsideCircle == Vertices.Count)
            {
                if (IsPointInsidePolygon(circle.Center)) return CirclePolygonPosition.CircleInsidePolygon;
                else return CirclePolygonPosition.Disjoint;
            }
            else if (verticesOutsideCircle == 0) return CirclePolygonPosition.PolygonInsideCircle;
            else if ((verticesOnCircle == 1) && (verticesInsideCircle == 0)) return CirclePolygonPosition.Disjoint;
            else return CirclePolygonPosition.Intersecting;
        }

        public PolygonPointPosition FindRelativePositionOfPoint(ICartesianPoint2D point)
        {
            
            int edgesPassingThroughPoint = 0;
            for (int i = 0; i < Vertices.Count; ++i)
            {
                var edge = new LineSegment2D(Vertices[i], Vertices[(i + 1) % Vertices.Count]); //TODO: Perhaps I should use the DirectedSegment2D, which is faster for the use
                if (edge.FindRelativePositionOfPoint(point) == LineSegment2D.SegmentPointPosition.PointOnSegment)
                {
                    ++edgesPassingThroughPoint;
                }
            }
            if (edgesPassingThroughPoint == 1) return PolygonPointPosition.OnEdge;
            else if (edgesPassingThroughPoint == 2) return PolygonPointPosition.OnVertex;
            else if (IsPointInsidePolygon(point)) return PolygonPointPosition.Inside;
            else 
            if (IsPointInsidePolygon(point)) return PolygonPointPosition.Inside;
            else return PolygonPointPosition.Outside;

            //TODO: IsPointInsidePolygon(point) is undefined the point is on the boundary, so the following doesn't always work.
            //if (IsPointInsidePolygon(point)) return PolygonPointPosition.Inside;
            //{
            //    int edgesPassingThroughPoint = 0;
            //    for (int i = 0; i < Vertices.Count; ++i)
            //    {
            //        var edge = new LineSegment2D(Vertices[i], Vertices[(i + 1) % Vertices.Count]); //TODO: Perhaps I should use the DirectedSegment2D, which is faster for the use
            //        if (edge.FindRelativePositionOfPoint(point) == LineSegment2D.SegmentPointPosition.PointOnSegment)
            //        {
            //            ++edgesPassingThroughPoint;
            //        }
            //    }
            //    if (edgesPassingThroughPoint == 0) return PolygonPointPosition.Outside;
            //    else if (edgesPassingThroughPoint == 1) return PolygonPointPosition.OnEdge;
            //    else if (edgesPassingThroughPoint == 2) return PolygonPointPosition.OnVertex;
            //    else throw new Exception("This should not have happened.");
            //}

        }

        /// <summary>
        /// Shortcut method to avoid redundant checks
        /// </summary>
        /// <param name="circle"></param>
        /// <returns></returns>
        public bool IntersectsWithCircle(Circle2D circle)
        {
            int verticesOutsideCircle = 0;
            int verticesInsideCircle = 0;
            foreach (ICartesianPoint2D vertex in Vertices)
            {
                CirclePointPosition vertexPosition = circle.FindRelativePositionOfPoint(vertex);
                if (vertexPosition == CirclePointPosition.Outside) ++verticesOutsideCircle;
                else if (vertexPosition == CirclePointPosition.Inside) ++verticesInsideCircle;
            }
            int verticesOnCircle = Vertices.Count - verticesOutsideCircle - verticesInsideCircle;

            if (verticesOutsideCircle >= 1)
            {
                if ((verticesInsideCircle >= 1) || (verticesOnCircle >= 2)) return true;
            }
            return false;
        }
    }
}
