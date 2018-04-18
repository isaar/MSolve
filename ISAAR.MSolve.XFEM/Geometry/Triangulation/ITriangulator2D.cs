using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Geometry.Triangulation
{
    interface ITriangulator2D
    {
        IReadOnlyList<Triangle2D> CreateMesh(IEnumerable<INaturalPoint2D> points);
        IReadOnlyList<Triangle2D> CreateMesh(IEnumerable<INaturalPoint2D> points, double maxTriangleArea);
    }
}
