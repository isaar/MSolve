using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TriangleNet;
using TriangleNet.Geometry;
using TriangleNet.Meshing;
using TriangleNet.Meshing.Algorithm;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;


namespace ISAAR.MSolve.XFEM.Geometry.Triangulation
{
    /// <summary>
    /// Wrapper class for 3rd party library incremental triangulator.
    /// </summary>
    class IncrementalTriangulator
    {
        private readonly List<Vertex> vertices;
        private readonly Incremental mesher;
        private readonly Configuration config;

        public IncrementalTriangulator(INaturalPoint2D vertex0, INaturalPoint2D vertex1, INaturalPoint2D vertex2)
        {
            this.vertices = new List<Vertex>();
            vertices.Add(new Vertex(vertex0.Xi, vertex0.Eta));
            vertices.Add(new Vertex(vertex1.Xi, vertex1.Eta));
            vertices.Add(new Vertex(vertex2.Xi, vertex2.Eta));

            this.mesher = new Incremental();
            this.config = new Configuration();
        }

        public void Insert(INaturalPoint2D vertex)
        {
            vertices.Add(new Vertex(vertex.Xi, vertex.Eta));
        }

        public IReadOnlyList<Triangle2D> Mesh()
        {
            var triangles = new List<Triangle2D>();
            IMesh mesh = mesher.Triangulate(vertices, config);
            foreach (ITriangle triangle in mesh.Triangles)
            {
                triangles.Add(new Triangle2D(triangle));
            }
            return triangles;
        }
    }
}
