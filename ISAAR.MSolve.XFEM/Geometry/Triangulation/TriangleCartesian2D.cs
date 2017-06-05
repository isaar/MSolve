using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Geometry.Triangulation
{
    // TODO: All triangles and shapes should be for cartesian system. They were previously for the Natural because I  
    // was following Nguyen's practice, which was actually incorrect.
    class TriangleCartesian2D
    {
        public IReadOnlyList<ICartesianPoint2D> Vertices { get; }

        public TriangleCartesian2D(TriangleNet.Geometry.ITriangle triangle)
        {
            var vertices = new CartesianPoint2D[3];
            for (int i = 0; i < 3; ++i)
            {
                vertices[i] = new CartesianPoint2D(triangle.GetVertex(i).X, triangle.GetVertex(i).Y);
            }
            this.Vertices = vertices;
        }

        public TriangleCartesian2D(CartesianPoint2D p0, CartesianPoint2D p1, CartesianPoint2D p2)
        {
            var vertices = new CartesianPoint2D[3];
            vertices[0] = p0;
            vertices[1] = p1;
            vertices[2] = p2;
            this.Vertices = vertices;
        }

        public override string ToString()
        {
            var builder = new StringBuilder("{ ");
            builder.Append(Vertices[0]).Append(" , ");
            builder.Append(Vertices[1]).Append(" , ");
            builder.Append(Vertices[2]);
            builder.Append(" }");
            return builder.ToString();
        }
    }
}
