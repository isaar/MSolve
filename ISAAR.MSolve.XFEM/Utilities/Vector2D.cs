using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Utilities
{
    class Vector2D
    {
        public double X { get; }
        public double Y { get; }
        public double Length { get; }

        public Vector2D(double xStart, double yStart, double xEnd, double yEnd)
        {
            X = xEnd - xStart;
            Y = yEnd - yStart;
            Length = Math.Sqrt(X * X + Y * Y);
        }

        public double DotProduct(Vector2D other)
        {
            return this.X * other.X + this.Y * other.Y;
        }

    }
}
