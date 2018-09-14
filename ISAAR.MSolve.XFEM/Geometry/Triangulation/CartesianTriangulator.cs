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
    // TODO: All triangulators should be for cartesian system. They were previously for the Natural because I was 
    // following Nguyen's practice, which was actually incorrect.
    class CartesianTriangulator
    {
        private readonly Dwyer mesher;
        private readonly Configuration config;

        public CartesianTriangulator()
        {
            this.mesher = new Dwyer();
            this.config = new Configuration();
        }

        public IReadOnlyList<TriangleCartesian2D> CreateMesh(IEnumerable<ICartesianPoint2D> points)
        {
            List<Vertex> vertices = new List<Vertex>();
            foreach (ICartesianPoint2D point in points)
            {
                vertices.Add(new Vertex(point.X, point.Y));
            }

            var triangles = new List<TriangleCartesian2D>();
            IMesh mesh = mesher.Triangulate(vertices, config);
            foreach (ITriangle triangle in mesh.Triangles)
            {
                triangles.Add(new TriangleCartesian2D(triangle));
            }
            return triangles;
        }
    }
}
