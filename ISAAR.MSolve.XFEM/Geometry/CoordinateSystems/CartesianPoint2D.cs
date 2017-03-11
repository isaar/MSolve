using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Geometry.CoordinateSystems
{
    class CartesianPoint2D: ICartesianPoint2D
    {
        public double X { get; }
        public double Y { get; }

        public CartesianPoint2D(double x, double y)
        {
            this.X = x;
            this.Y = y;
        }

        public override string ToString()
        {
            return "(x , y) = (" + X + " , " + Y + ")";
        }
    }
}
