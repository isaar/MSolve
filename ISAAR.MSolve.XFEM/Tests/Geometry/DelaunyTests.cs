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
using ISAAR.MSolve.XFEM.Geometry.Triangulation;

namespace ISAAR.MSolve.XFEM.Tests.Geometry
{
    class DelaunyTests
    {
        private static readonly NaturalPoint2D[] nodes =
        {
            new NaturalPoint2D(-1.0, -1.0),
            new NaturalPoint2D(1.0, -1.0),
            new NaturalPoint2D(1.0, 1.0),
            new NaturalPoint2D(-1.0, 1.0)
        };

        private static readonly NaturalPoint2D[] intersectionPoints1 =
        {
            new NaturalPoint2D(-1.0, 0.0),
            new NaturalPoint2D(1.0, 0.0)
        };

        private static readonly NaturalPoint2D[] intersectionPoints2 =
        {
            new NaturalPoint2D(-1.0, 0.0),
            new NaturalPoint2D(-1.0/3, 0.0)
        };

        private static void TestMesh(NaturalPoint2D[] intersectionPoints)
        {
            var mesher = new IncrementalTriangulator();
            var vertices = new List<INaturalPoint2D>();
            vertices.AddRange(nodes);
            vertices.AddRange(intersectionPoints);
            var triangles = mesher.CreateMesh(vertices);
            Utilities.Triangles.PrintMesh(triangles);
        }

        public static void Main()
        {
            TestMesh(intersectionPoints1);
            TestMesh(intersectionPoints2);
        }
    }
}
