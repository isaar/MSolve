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

        // The local coordinates of the tip.
        private readonly double x1Tip;
        private readonly double x2Tip;

        public TipCoordinateSystem(ICartesianPoint2D tipCoordinates, double orientation)
        {
            cosa = Math.Cos(orientation);
            sina = Math.Sin(orientation);
            x1Tip = cosa * tipCoordinates.X + sina * tipCoordinates.Y;
            x2Tip = -sina * tipCoordinates.X + cosa * tipCoordinates.Y;
        }

        public PolarPoint2D TransformCartesianGlobalToPolarLocal(ICartesianPoint2D cartesianGlobalPoint)
        {
            double globX = cartesianGlobalPoint.X;
            double globY = cartesianGlobalPoint.Y;
            double x1 = cosa * globX + sina * globY - x1Tip;
            double x2 = -sina * globX + cosa * globY - x2Tip;
            double r = Math.Sqrt(x1 * x1 + x2 * x2);
            double theta = Math.Atan2(x2, x1);
            #region Testing
            // To compare with C++ code, where ATAN2 is [-π, π) uncomment the next following
            //if (Math.Abs(x2) < 1e-8) theta = -theta;
            #endregion
            return new PolarPoint2D(r, theta);
        }

        public TipJacobian CalculateJacobianAt(PolarPoint2D polarCoordinates)
        {
            return new TipJacobian(cosa, sina, polarCoordinates.R, polarCoordinates.Theta);
        }
    }
}
