using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Geometry.Shapes
{
    class ConvexPolygon2D
    {
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

        public IReadOnlyList<ICartesianPoint2D> Vertices { get; }
        public IReadOnlyList<LineSegment2D> Edges { get; }

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
    }
}
