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
    class IncrementalTriangulator: ITriangulator2D
    {
        private readonly Incremental mesher;
        private readonly Configuration config;

        public IncrementalTriangulator()
        {
            this.mesher = new Incremental();
            this.config = new Configuration();
        }

        public IReadOnlyList<Triangle2D> CreateMesh(IEnumerable<INaturalPoint2D> points)
        {
            List<Vertex> vertices = new List<Vertex>();
            foreach (INaturalPoint2D point in points)
            {
                vertices.Add(new Vertex(point.Xi, point.Eta));
            }

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
