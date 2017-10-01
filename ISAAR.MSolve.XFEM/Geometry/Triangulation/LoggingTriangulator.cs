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
    class LoggingTriangulator : ITriangulator2D
    {
        private readonly Incremental mesher;
        private readonly Configuration config;
        private readonly int elementID;

        public LoggingTriangulator(int elementID)
        {
            this.mesher = new Incremental();
            this.config = new Configuration();
            this.elementID = elementID;
        }

        public IReadOnlyList<Triangle2D> CreateMesh(IEnumerable<INaturalPoint2D> points)
        {
            List<Vertex> vertices = new List<Vertex>();
            foreach (INaturalPoint2D point in points)
            {
                vertices.Add(new Vertex(point.Xi, point.Eta));
            }

            IMesh mesh = mesher.Triangulate(vertices, config);
            var triangles = new List<Triangle2D>();
            foreach (ITriangle triangle in mesh.Triangles)
            {
                triangles.Add(new Triangle2D(triangle));
            }

            Console.WriteLine("********* Interjection from triangulator of element " + elementID + " - START ************");
            Console.WriteLine("Generated Mesh:");
            Utilities.Triangles.PrintMesh(triangles);
            Console.WriteLine("********* Interjection from triangulator of element " + elementID + " - END ************\n");

            return triangles;
        }

        public IReadOnlyList<Triangle2D> CreateMesh(IEnumerable<INaturalPoint2D> points, double maxTriangleArea)
        {
            var quality = new QualityOptions();
            quality.MaximumArea = maxTriangleArea;

            List<Vertex> vertices = new List<Vertex>();
            foreach (INaturalPoint2D point in points)
            {
                vertices.Add(new Vertex(point.Xi, point.Eta));
            }

            IMesh mesh = mesher.Triangulate(vertices, config);
            mesh.Refine(quality, true);

            var triangles = new List<Triangle2D>();
            foreach (ITriangle triangle in mesh.Triangles)
            {
                triangles.Add(new Triangle2D(triangle));
            }

            Console.WriteLine("********* Interjection from triangulator of element " + elementID + " - START ************");
            Console.WriteLine("Generated Mesh:");
            Utilities.Triangles.PrintMesh(triangles);
            Console.WriteLine("********* Interjection from triangulator of element " + elementID + " - END ************\n");

            return triangles;
        }
    }
}
