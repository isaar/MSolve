using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry;

namespace ISAAR.MSolve.XFEM.Interpolation
{
    class ShapeFunctionValues2D
    {
        private readonly double[] values;
        private readonly IReadOnlyList<Node2D> nodes;

        public ShapeFunctionValues2D(double[] shapeFunctionValues, IReadOnlyList<Node2D> nodes)
        {
            this.values = shapeFunctionValues;
            this.nodes = nodes;
        }

        public double this[int nodeIndex] { get { return values[nodeIndex]; } }

        public IPoint2D TransformNaturalToCartesian(IPoint2D naturalCoordinates)
        {
            double x = 0, y = 0;
            for (int nodeIdx = 0; nodeIdx < values.Length; ++nodeIdx)
            {
                x += values[nodeIdx] * nodes[nodeIdx].X;
                y += values[nodeIdx] * nodes[nodeIdx].Y;
            }
            return new Point2D(x, y);
        }
    }
}
