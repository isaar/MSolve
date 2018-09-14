using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Providers;

namespace ISAAR.MSolve.XFEM.Tests.Geometry
{
    class RectilinearTest
    {
        public static void Main()
        {
            double[,] xSizes = new double[,] { { 0, 1 }, { 4, 0.5 }, { 7, 0.1 }, { 9, 0.5 }, { 13, 1 } };
            var generator = new RectilinearMeshGenerator(xSizes, xSizes);
            var meshEntities = generator.CreateMesh();
            foreach (var node in meshEntities.Item1)
            {
                Console.WriteLine(node);
            }
        }
    }
}
