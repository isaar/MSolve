using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry;

namespace ISAAR.MSolve.XFEM.Integration.Points
{
    class GaussPoint2D: IPoint2D
    {
        public double X { get; }
        public double Y { get; }
        public double Weight { get; }

        public GaussPoint2D(double x, double y, double weight)
        {
            this.X = x;
            this.Y = y;
            this.Weight = weight;
        }
    }
}
