using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.Tests.Geometry
{
    class PointInPolygon
    {
        private static readonly ConvexPolygon2D rectangle = ConvexPolygon2D.CreateUnsafe(new CartesianPoint2D[]
            {
                new CartesianPoint2D(-4.0, -2.5),
                new CartesianPoint2D(3.0, -2.5),
                new CartesianPoint2D(3.0, 1.0),
                new CartesianPoint2D(-4.0, 1.0),
            });

        private static void PointInsidePolygon()
        {
            var point = new CartesianPoint2D(1.75, -0.25);
            if (rectangle.IsPointInsidePolygon(point)) Console.WriteLine("Point is inside polygon, as expected");
            else Console.WriteLine("Error: point is outside polygon");
        }

        private static void PointOutsidePolygon()
        {
            var point = new CartesianPoint2D(10.0, -0.25);
            if (!rectangle.IsPointInsidePolygon(point)) Console.WriteLine("Point is outside polygon, as expected");
            else Console.WriteLine("Error: point is inside polygon");
        }

        private static void PointOnEdgeOfPolygon()
        {
            var point = new CartesianPoint2D(1.0, 1.0);
            if (rectangle.IsPointInsidePolygon(point)) Console.WriteLine("Point is inside polygon");
            else Console.WriteLine("Point is outside polygon");
        }

        private static void PointBeforeEdgeOfPolygon()
        {
            var point = new CartesianPoint2D(-10.0, 1.0);
            if (rectangle.IsPointInsidePolygon(point)) Console.WriteLine("Point is inside polygon");
            else Console.WriteLine("Point is outside polygon");
        }

        private static void PointAfterEdgeOfPolygon()
        {
            var point = new CartesianPoint2D(10.0, 1.0);
            if (rectangle.IsPointInsidePolygon(point)) Console.WriteLine("Point is inside polygon");
            else Console.WriteLine("Point is outside polygon");
        }

        private static void PointOnVertexOfPolygon()
        {
            for (int i = 0; i < rectangle.Vertices.Count; ++i)
            {
                if (rectangle.IsPointInsidePolygon(rectangle.Vertices[i]))
                    Console.WriteLine("Vertex " + i + " is inside polygon");
                else Console.WriteLine("Vertex " + i + " is outside polygon");
            }
            
        }

        public static void Main()
        {
            //PointInsidePolygon();
            //PointOutsidePolygon();
            PointOnEdgeOfPolygon();
            PointBeforeEdgeOfPolygon();
            PointAfterEdgeOfPolygon();
            PointOnVertexOfPolygon();
        }
    }
}
