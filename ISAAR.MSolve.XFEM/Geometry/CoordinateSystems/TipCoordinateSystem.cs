using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.XFEM.Geometry.CoordinateSystems
{
    // Isn't the orientation the same as the polar coordinate theta?
    class TipCoordinateSystem
    {
        // Let a be the counter-clockwise angle from global x axis to local x1 axis
        public readonly double cosa;
        public readonly double sina;

        // The local coordinates of the global origin
        private readonly double x1Origin;
        private readonly double x2Origin;

        public TipCoordinateSystem(ICartesianPoint2D tipCoordinates, double orientation)
        {
            cosa = Math.Cos(orientation);
            sina = Math.Sin(orientation);
            x1Origin = - cosa * tipCoordinates.X - sina * tipCoordinates.Y;
            x2Origin = sina * tipCoordinates.X - cosa * tipCoordinates.Y;
        }

        public PolarPoint2D TransformCartesianGlobalToPolarLocal(ICartesianPoint2D cartesianGlobalPoint)
        {
            double globX = cartesianGlobalPoint.X;
            double globY = cartesianGlobalPoint.Y;
            double x1 = cosa * globX + sina * globY + x1Origin;
            double x2 = -sina * globX + cosa * globY + x2Origin;
            double r = Math.Sqrt(x1 * x1 + x2 * x2);
            double theta = Math.Atan2(x2, x1);
            return new PolarPoint2D(r, theta);
        }

        public TipJacobian CalculateJacobianAt(PolarPoint2D polarCoordinates)
        {
            return new TipJacobian(cosa, sina, polarCoordinates.R, polarCoordinates.Theta);
        }
    }
}
