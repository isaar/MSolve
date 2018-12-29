using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;

namespace ISAAR.MSolve.XFEM.Utilities
{
    static class Triangles
    {
        public static void PrintMesh(IReadOnlyList<Triangle2D> triangles)
        {
            for (int i = 0; i < triangles.Count; ++i)
            {
                Console.Write("Triangle " + i + ": ");
                Console.WriteLine(triangles[i]);
            }
        }
    }
}
