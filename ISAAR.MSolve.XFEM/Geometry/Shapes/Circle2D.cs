using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.XFEM.Geometry.Shapes
{
    public enum CirclePointPosition
    {
        Inside, On, Outside
    }

    class Circle2D
    {
        public CartesianPoint2D Center { get; }
        public double Radius { get; }

        public Circle2D(CartesianPoint2D center, double radius)
        {
            this.Center = new CartesianPoint2D(center.X, center.Y); // Copy it for extra safety
            this.Radius = radius;
        }

        public CirclePointPosition FindRelativePositionOfPoint(CartesianPoint2D point)
        {
            double distanceX = point.X - Center.X;
            double distanceY = point.Y - Center.Y;
            double distance = Math.Sqrt(distanceX * distanceX + distanceY * distanceY);

            if (distance > Radius) return CirclePointPosition.Outside;
            else if (distance < Radius) return CirclePointPosition.Inside;
            else return CirclePointPosition.On;
        }
    }
}
