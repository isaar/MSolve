using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Geometry
{
    class Point2D : IPoint2D
    {
        public double X { get; }
        public double Y { get; }

        public Point2D(double x, double y)
        {
            this.X = x;
            this.Y = y;
        }
    }
}
