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
    class MirrorNaturalXTriangulator : ITriangulator2D
    {
        private readonly Incremental mesher;
        private readonly Configuration config;
        private readonly int elementID;

        public MirrorNaturalXTriangulator(int elementID)
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

            var triangles = new List<Triangle2D>();
            IMesh mesh = mesher.Triangulate(vertices, config);
            foreach (ITriangle triangle in mesh.Triangles)
            {
                double xi0 = triangle.GetVertex(0).X;
                double eta0 = -triangle.GetVertex(0).Y;
                double xi1 = triangle.GetVertex(1).X;
                double eta1 = -triangle.GetVertex(1).Y;
                double xi2 = triangle.GetVertex(2).X;
                double eta2 = -triangle.GetVertex(2).Y;
                triangles.Add(new Triangle2D(
                    new NaturalPoint2D(xi0, eta0), new NaturalPoint2D(xi1, eta1), new NaturalPoint2D(xi2, eta2)));
            }

            //Console.WriteLine("********* Interjection from triangulator of element " + elementID + " - START ************");
            //Console.WriteLine("Generated Mesh:");
            //Utilities.Triangles.PrintMesh(triangles);
            //Console.WriteLine();
            //Console.WriteLine("********* Interjection from triangulator of element " + elementID + " - END ************\n");

            return triangles;
        }

        public IReadOnlyList<Triangle2D> CreateMesh(IEnumerable<INaturalPoint2D> points, double maxTriangleArea)
        {
            throw new NotImplementedException();
        }
    }
}
