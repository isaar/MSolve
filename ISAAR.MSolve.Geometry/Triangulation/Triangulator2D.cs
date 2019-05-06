using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TriangleNet;
using TriangleNet.Geometry;
using TriangleNet.Meshing;
using TriangleNet.Meshing.Algorithm;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.Geometry.Triangulation
{
    // TODO: All triangulators should be for cartesian system. They were previously for the Natural because I was 
    // following Nguyen's practice, which was actually incorrect.
    public class Triangulator2D<TVertex> : ITriangulator2D<TVertex> where TVertex : IPoint
    {
        private readonly CreateVertex2D<TVertex> createVertex;
        private readonly Configuration config;
        private readonly TriangleNet.Meshing.ITriangulator mesher;

        public Triangulator2D(CreateVertex2D<TVertex> createVertex) : this(new Dwyer(), createVertex) {}

        /// <summary>
        /// 
        /// </summary>
        /// <param name="mesher">The Triangle.NET algorithm: eg. Dwyer, Incremental</param>
        /// <param name="createVertex"></param>
        public Triangulator2D(TriangleNet.Meshing.ITriangulator mesher, CreateVertex2D<TVertex> createVertex)
        {
            this.createVertex = createVertex;
            this.config = new Configuration();
            this.mesher = mesher;
        }

        public IReadOnlyList<Triangle2D<TVertex>> CreateMesh(IEnumerable<TVertex> points)
        {
            List<Vertex> vertices = new List<Vertex>();
            foreach (TVertex point in points) vertices.Add(new Vertex(point.X1, point.X2));

            var triangles = new List<Triangle2D<TVertex>>();
            IMesh mesh = mesher.Triangulate(vertices, config);
            foreach (ITriangle triangle in mesh.Triangles) triangles.Add(new Triangle2D<TVertex>(triangle, createVertex));
            return triangles;
        }

        public IReadOnlyList<Triangle2D<TVertex>> CreateMesh(IEnumerable<TVertex> points, double maxTriangleArea)
        {
            // TODO: Refinement gives incorrect results. Probably the refined mesh doesn't conform to crack.
            throw new NotImplementedException("Refinement gives incorrect results. " +
                "Probably the refined mesh doesn't conform to crack.");

            //var quality = new QualityOptions();
            //quality.MaximumArea = maxTriangleArea;

            //var vertices = new List<TriangleNet.Geometry.Vertex>();
            //foreach (TVertex point in points) vertices.Add(new TriangleNet.Geometry.Vertex(point.X1, point.X2));

            //IMesh mesh = mesher.Triangulate(vertices, config);
            //mesh.Refine(quality, true);

            //var triangles = new List<Triangle2D<TVertex>>();
            //foreach (ITriangle triangle in mesh.Triangles) triangles.Add(new Triangle2D<TVertex>(triangle, createVertex));
            //return triangles;
        }
    }
}
