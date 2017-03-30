using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Utilities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Enrichments.Items.CrackTip
{
    // Perhaps I should not expose this class, but use it privately in TipCoordinateSystem and have batch methods when 
    // a lot of calls need to calculate the jacobians at the same point.
    class TipJacobians
    {
        private readonly TipCoordinateSystem tipSystem;
        private readonly double[,] invTransJacobianPolarToLocal;
        private readonly double[,] invTransJacobianPolarToGlobal;

        public TipJacobians(TipCoordinateSystem tipSystem, PolarPoint2D polarCoordinates)
        {
            this.tipSystem = tipSystem;

            double r = polarCoordinates.R;
            double cosTheta = Math.Cos(polarCoordinates.Theta);
            double sinTheta = Math.Sin(polarCoordinates.Theta);
            invTransJacobianPolarToLocal = new double[,] { { cosTheta, -sinTheta / r }, { sinTheta, cosTheta / r } };

            invTransJacobianPolarToGlobal = 
                tipSystem.TransposeRotationMatrixGlobalToLocal.MultiplyRight(invTransJacobianPolarToLocal);
        }

        public Tuple<double, double> TransformScalarFieldDerivativesLocalPolarToLocalCartesian(
            double gradR, double gradTheta)
        {
            double gradX1 = invTransJacobianPolarToLocal[0, 0] * gradR
                + invTransJacobianPolarToLocal[0, 1] * gradTheta;
            double gradX2 = invTransJacobianPolarToLocal[1, 0] * gradR
                +  invTransJacobianPolarToLocal[1, 1] * gradTheta;
           return new Tuple<double, double>(gradX1, gradX2);
        }

        public Tuple<double, double> TransformScalarFieldDerivativesLocalPolarToGlobalCartesian(
            double derivativeR, double derivativeTheta)
        {
            double gradX1 = invTransJacobianPolarToGlobal[0, 0] * derivativeR
                + invTransJacobianPolarToGlobal[0, 1] * derivativeTheta;
            double gradX2 = invTransJacobianPolarToGlobal[1, 0] * derivativeR
                + invTransJacobianPolarToGlobal[1, 1] * derivativeTheta;
            return new Tuple<double, double>(gradX1, gradX2);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="derivatives">Row i is the gradient of the ith component of the vector field, thus: 
        ///     derivatives = [Fr,r Fr,theta; Ftheta,r Ftheta,theta],
        ///     where Fi,j is the derivative of component i w.r.t. coordinate j</param>
        /// <returns></returns>
        public double[,] TransformVectorFieldDerivativesLocalPolarToLocalCartesian(double[,] derivatives)
        {
            throw new NotImplementedException("Should return [derivatives]*inv([JacobianPolarToLocal])");
        }
    }
}
