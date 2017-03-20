using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Geometry.Triangulation
{
    // 
    /// <summary>
    /// Wrapper class for 3rd party library triangles. For now I only need 2D triangles in natural system. 
    /// TODO: make them generic.
    /// </summary>
    class Triangle2D
    {
        public IReadOnlyList<INaturalPoint2D> Vertices { get; }

        public Triangle2D(TriangleNet.Geometry.ITriangle triangle)
        {
            var vertices = new NaturalPoint2D[3];
            for (int i = 0; i < 3; ++i)
            {
                vertices[i] = new NaturalPoint2D(triangle.GetVertex(i).X, triangle.GetVertex(i).Y);
            }
            this.Vertices = vertices;
        }

        public Triangle2D(NaturalPoint2D p0, NaturalPoint2D p1, NaturalPoint2D p2)
        {
            var vertices = new NaturalPoint2D[3];
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
