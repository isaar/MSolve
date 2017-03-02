using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.XFEM.Geometry.CoordinateSystems
{
    class TipJacobian
    {
        // TODO: Purge these after testing
        private readonly double cosa, sina;
        private readonly Matrix2D<double> jacobian;

        public TipJacobian(double cosa, double sina, double r, double theta)
        {
            this.cosa = cosa;
            this.sina = sina;

            double cosTheta = Math.Cos(theta);
            double sinTheta = Math.Sin(theta);

            jacobian = new Matrix2D<double>(2, 2);
            jacobian[0, 0] = cosa * cosTheta - sina * sinTheta;
            jacobian[0, 1] = - cosa * sinTheta / r - sina * cosTheta / r;
            jacobian[1, 0] = sina * cosTheta + cosa * sinTheta;
            jacobian[1, 1] = - sina * sinTheta / r + cosa * cosTheta / r;
        }

        public Tuple<double, double> TransformLocalPolarDerivativesToGlobalCartesian(
            double derivativeR, double derivativeTheta)
        {
            double derivativeX = derivativeR * jacobian[0, 0] + derivativeTheta * jacobian[0, 1];
            double derivativeY = derivativeR * jacobian[1, 0] + derivativeTheta * jacobian[1, 1];
            return new Tuple<double, double>(derivativeX, derivativeY);
        }

        public Tuple<double, double> TransformLocalCartesianDerivativesToGlobalCartesian(
            double derivativeX1, double derivativeX2)
        {
            double derivativeX = derivativeX1 * cosa - derivativeX2 * sina;
            double derivativeY = derivativeX1 * sina + derivativeX2 * cosa;
            return new Tuple<double, double>(derivativeX, derivativeY);
        }
    }
}
