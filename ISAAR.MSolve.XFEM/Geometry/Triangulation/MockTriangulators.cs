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
    class MockTriangulators
    {
        public static readonly ITriangulator2D ForElement4 = new Element4Triangulator();
        public static readonly ITriangulator2D ForElement5 = new Element5Triangulator();

        private class Element4Triangulator: ITriangulator2D
        {
            public IReadOnlyList<Triangle2D> CreateMesh(IEnumerable<INaturalPoint2D> points)
            {
                var triangles = new Triangle2D[4];
                triangles[0] = new Triangle2D(
                    new NaturalPoint2D(1, -0.15), new NaturalPoint2D(1, 1), new NaturalPoint2D(-1, 1));
                triangles[1] = new Triangle2D(
                    new NaturalPoint2D(-1, -0.599999999999999), new NaturalPoint2D(1, -0.15), new NaturalPoint2D(-1, 1));
                triangles[2] = new Triangle2D(
                    new NaturalPoint2D(-1, -0.599999999999999), new NaturalPoint2D(-1, -1), new NaturalPoint2D(1, -1));
                triangles[3] = new Triangle2D(
                    new NaturalPoint2D(-1, -0.599999999999999), new NaturalPoint2D(1, -1), new NaturalPoint2D(1, -0.15));
                return triangles;
            }

            public IReadOnlyList<Triangle2D> CreateMesh(IEnumerable<INaturalPoint2D> points, double maxTriangleArea)
            {
                throw new NotImplementedException();
            }
        }

        private class Element5Triangulator : ITriangulator2D
        {
            public IReadOnlyList<Triangle2D> CreateMesh(IEnumerable<INaturalPoint2D> points)
            {
                var triangles = new Triangle2D[5];
                triangles[0] = new Triangle2D(
                    new NaturalPoint2D(-0.333333333333333, 5.55111512312578e-16), new NaturalPoint2D(-1, -1), new NaturalPoint2D(1, -1));
                triangles[1] = new Triangle2D(
                    new NaturalPoint2D(-0.333333333333333, 5.55111512312578e-16), new NaturalPoint2D(-1, -0.15), new NaturalPoint2D(-1, -1));
                triangles[2] = new Triangle2D(
                    new NaturalPoint2D(-0.333333333333333, 5.55111512312578e-16), new NaturalPoint2D(1, -1), new NaturalPoint2D(1, 1));
                triangles[3] = new Triangle2D(
                    new NaturalPoint2D(-0.333333333333333, 5.55111512312578e-16), new NaturalPoint2D(1, 1), new NaturalPoint2D(-1, 1));
                triangles[4] = new Triangle2D(
                    new NaturalPoint2D(-0.333333333333333, 5.55111512312578e-16), new NaturalPoint2D(-1, 1), new NaturalPoint2D(-1, -0.15));
                return triangles;
            }

            public IReadOnlyList<Triangle2D> CreateMesh(IEnumerable<INaturalPoint2D> points, double maxTriangleArea)
            {
                throw new NotImplementedException();
            }
        }

    }
}