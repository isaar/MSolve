using System.Collections.Generic;
using ISAAR.MSolve.Geometry.Commons;

namespace ISAAR.MSolve.Geometry.Coordinates
{
    public class Point2DComparerXMajor: IComparer<CartesianPoint>
    {
        private readonly ValueComparer valueComparer;

        public Point2DComparerXMajor(double tolerance = 1e-6)
        {
            this.valueComparer = new ValueComparer(tolerance);
        }

        public int Compare(CartesianPoint point1, CartesianPoint point2)
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
