using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Interpolation.InverseMappings
{
    class InverseTri3Mapping: IInverseMapping2D
    {
        private readonly double x1, x2, x3, y1, y2, y3;
        private readonly double det;

        public InverseTri3Mapping(IReadOnlyList<Node2D> nodes)
        {
            x1 = nodes[0].X;
            x2 = nodes[1].X;
            x3 = nodes[2].X;
            y1 = nodes[0].Y;
            y2 = nodes[1].Y;
            y3 = nodes[2].Y;
            det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
        }

        public INaturalPoint2D TransformCartesianToNatural(ICartesianPoint2D point)
        {
            double detXi = (point.X - x1) * (y3 - y1) - (x3 - x1) * (point.Y - y1);
            double detEta = (x2 - x1) * (point.Y - y1) - (point.X - x1) * (y2 - y1);
            return new NaturalPoint2D(detXi / det, detEta / det);
        }
    }
}
