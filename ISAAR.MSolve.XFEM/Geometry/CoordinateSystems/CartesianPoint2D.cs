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

        public CartesianPoint2D(Vector2 cartesianCoordinates)
        {
            this.X = cartesianCoordinates[0];
            this.Y = cartesianCoordinates[1];
        }

        public Vector2 Coordinates { get { return Vector2.Create( X, Y ); } }

        public override string ToString()
        {
            return "(x , y) = (" + X + " , " + Y + ")";
        }
    }
}
