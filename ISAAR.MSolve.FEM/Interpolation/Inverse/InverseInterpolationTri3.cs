using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.FEM.Interpolation.Inverse
{
    /// <summary>
    /// Inverse mapping of the isoparametric interpolation of a triangular finite element with 3 nodes. Since the original 
    /// mapping is linear, there are analytic formulas.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class InverseInterpolationTri3 : IInverseInterpolation2D
    {
        private readonly double x1, x2, x3, y1, y2, y3;
        private readonly double det;

        public InverseInterpolationTri3(IReadOnlyList<Node2D> nodes)
        {
            x1 = nodes[0].X;
            x2 = nodes[1].X;
            x3 = nodes[2].X;
            y1 = nodes[0].Y;
            y2 = nodes[1].Y;
            y3 = nodes[2].Y;
            det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
        }

        public NaturalPoint2D TransformPointCartesianToNatural(CartesianPoint2D point)
        {
            double detXi = (point.X - x1) * (y3 - y1) - (x3 - x1) * (point.Y - y1);
            double detEta = (x2 - x1) * (point.Y - y1) - (point.X - x1) * (y2 - y1);
            return new NaturalPoint2D(detXi / det, detEta / det);
        }
    }
}
