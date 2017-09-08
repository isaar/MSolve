using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Geometry.CoordinateSystems
{
    class Point2DComparerXMajor: IComparer<ICartesianPoint2D>
    {
        public int Compare(ICartesianPoint2D point1, ICartesianPoint2D point2)
        {
            if (point1.X < point2.X) return -1;
            else if (point1.X > point2.X) return 1;
            else // same X
            {
                if (point1.Y < point2.Y) return -1;
                else if (point1.Y > point2.Y) return 1;
                else return 0; // same point
            }
        }
    }
}
