using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Commons;

namespace ISAAR.MSolve.Geometry.Coordinates
{
    public class Point2DComparerXMajor: IComparer<CartesianPoint2D>
    {
        private readonly ValueComparer valueComparer;

        public Point2DComparerXMajor(double tolerance = 1e-6)
        {
            this.valueComparer = new ValueComparer(tolerance);
        }

        public int Compare(CartesianPoint2D point1, CartesianPoint2D point2)
        {
            if (valueComparer.AreEqual(point1.X, point2.X))
            {
                if (valueComparer.AreEqual(point1.Y, point2.Y)) return 0;
                else if (point1.Y < point2.Y) return -1;
                else return 1;
            }
            else if (point1.X < point2.X) return -1;
            else return 1;
        }
    }
}
