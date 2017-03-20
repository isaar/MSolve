using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.XFEM.Geometry.CoordinateSystems
{
    class TipJacobian
    {
        private readonly Matrix2D jacobian;

        public TipJacobian(double cosa, double sina, double r, double theta)
        {
            double cosTheta = Math.Cos(theta);
            double sinTheta = Math.Sin(theta);

            jacobian = new Matrix2D(2, 2);
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
    }
}
