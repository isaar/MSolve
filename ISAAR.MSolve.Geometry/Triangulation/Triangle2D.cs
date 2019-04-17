using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.Geometry.Triangulation
{
    /// <summary>
    /// Initializes a new instance of some class that implements <see cref="IPoint"/>.
    /// </summary>
    /// <param name="x">The x coordinate of the node.</param>
    /// <param name="y">The y coordinate of the node.</param>
    public delegate TVertex CreateVertex2D<TVertex>(double x, double y) where TVertex : IPoint;

    /// <summary>
    /// Wrapper class for 3rd party library triangles. For now I only need 2D triangles in natural system. 
    /// TODO: make them generic.
    /// </summary>
    public class Triangle2D<TVertex> where TVertex: IPoint
    {
        public IReadOnlyList<TVertex> Vertices { get; }

        public Triangle2D(TriangleNet.Geometry.ITriangle triangle, CreateVertex2D<TVertex> createVertex)
        {
            var vertices = new TVertex[3];
            for (int i = 0; i < 3; ++i) vertices[i] = createVertex(triangle.GetVertex(i).X, triangle.GetVertex(i).Y);
            this.Vertices = vertices;
        }

        public Triangle2D(TVertex p0, TVertex p1, TVertex p2)
        {
            var vertices = new TVertex[3];
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
