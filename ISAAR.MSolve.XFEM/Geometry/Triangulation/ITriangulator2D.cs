using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.XFEM.Geometry.Triangulation
{
    interface ITriangulator2D
    {
        IReadOnlyList<Triangle2D> CreateMesh(IEnumerable<NaturalPoint2D> points);
        IReadOnlyList<Triangle2D> CreateMesh(IEnumerable<NaturalPoint2D> points, double maxTriangleArea);
    }
}
