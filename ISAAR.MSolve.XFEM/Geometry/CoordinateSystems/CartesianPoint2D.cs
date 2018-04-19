using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;

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

        public CartesianPoint2D(Vector coordinates)
        {
            this.X = coordinates[0];
            this.Y = coordinates[1];
        }

        public Vector Coordinates { get { return Vector.CreateFromArray(new double[] { X, Y }); } }

        public override string ToString()
        {
            return "(x , y) = (" + X + " , " + Y + ")";
        }
    }
}
