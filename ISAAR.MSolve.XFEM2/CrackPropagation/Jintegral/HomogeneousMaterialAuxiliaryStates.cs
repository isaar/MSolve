using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Tensors;

namespace ISAAR.MSolve.XFEM.CrackPropagation.Jintegral
{
    // TODO 1: An IHomogeneousMaterialField would be better than the restrictive HomogeneousElasticMaterial2D
    // TODO 2: Enforce that all elements of the (J-integral) domain have the identical material properties to the ones 
    //      passed to this class.
    class HomogeneousMaterialAuxiliaryStates : IAuxiliaryStates
    {
        private readonly HomogeneousElasticMaterial2D material;

        /// <summary>
        /// The material properties (E, v, E*, v*) must be the same across all elements. The user assumes 
        /// responsibility for passing a <see cref="HomogeneousElasticMaterial2D"/> that has the same properties as 
        /// the materials of all other elements of the integration domain.
        /// </summary>
        /// <param name="globalMaterial">The material properties which must be identical for all elements and this 
        ///     class</param>
        public HomogeneousMaterialAuxiliaryStates(HomogeneousElasticMaterial2D globalMaterial)
        {
            this.material = globalMaterial;
        }

        public AuxiliaryStatesTensors ComputeTensorsAt(ICartesianPoint2D globalIntegrationPoint, 
            TipCoordinateSystem tipCoordinateSystem)
        {
            // Common calculations
            PolarPoint2D polarCoordinates =
                tipCoordinateSystem.TransformPointGlobalCartesianToLocalPolar(globalIntegrationPoint);
            var commonValues = new CommonValues(polarCoordinates.R, polarCoordinates.Theta);

            // Displacement field derivatives
            Matrix2by2 displacementGradientMode1, displacementGradientMode2;
            TipJacobians polarJacobians = tipCoordinateSystem.CalculateJacobiansAt(polarCoordinates);
            ComputeDisplacementDerivatives(polarJacobians, commonValues, 
                out displacementGradientMode1, out displacementGradientMode2);

            // Strains
            Tensor2D strainTensorMode1 = ComputeStrainTensor(displacementGradientMode1);
            Tensor2D strainTensorMode2 = ComputeStrainTensor(displacementGradientMode2);

            // Stresses
            Tensor2D stressTensorMode1, stressTensorMode2;
            ComputeStressTensors(commonValues, out stressTensorMode1, out stressTensorMode2);

            return new AuxiliaryStatesTensors(displacementGradientMode1, displacementGradientMode2,
                strainTensorMode1, strainTensorMode2, stressTensorMode1, stressTensorMode2);
        }

        private void ComputeDisplacementDerivatives(TipJacobians polarJacobians, CommonValues val, 
            out Matrix2by2 displacementGradientMode1, out Matrix2by2 displacementGradientMode2)
        {
            // Temporary values and derivatives of the differentiated quantities. See documentation for their derivation.
            double E = material.HomogeneousYoungModulus;
            double v = material.HomogeneousPoissonRatio;
            double vEq = material.HomogeneousEquivalentPoissonRatio;
            double k = (3.0 - vEq) / (1.0 + vEq);
            double a = (1.0 + v) / (E * Math.Sqrt(2.0 * Math.PI));
            double b = val.sqrtR;
            double b_r = 0.5 / val.sqrtR;

            // Mode 1
            {
                // Temporary values that differ between the 2 modes
                double c1 = val.cosThetaOver2 * (k - val.cosTheta);
                double c2 = val.sinThetaOver2 * (k - val.cosTheta);
                double c1_theta = -0.5 * c2 + val.cosThetaOver2 * val.sinTheta;
                double c2_theta = 0.5 * c1 + val.sinThetaOver2 * val.sinTheta;

                // The vector field derivatives w.r.t. to the local polar coordinates. 
                // The vector components refer to the local cartesian system though.
                var polarGradient = Matrix2by2.CreateFromArray(
                    new double[,] { { a * b_r * c1, a * b * c1_theta }, { a * b_r * c2, a * b * c2_theta } });
                // The vector field derivatives w.r.t. to the local cartesian coordinates.
                displacementGradientMode1 = 
                    polarJacobians.TransformVectorFieldDerivativesLocalPolarToLocalCartesian(polarGradient);
            }

            // Mode 2
            {
                double paren1 = 2.0 + k + val.cosTheta;
                double paren2 = 2.0 - k - val.cosTheta;
                double c1 = val.sinThetaOver2 * paren1;
                double c2 = val.cosThetaOver2 * paren2;
                double c1_theta = 0.5 * val.cosThetaOver2 * paren1 - val.sinThetaOver2 * val.sinTheta;
                double c2_theta = -0.5 * val.sinThetaOver2 * paren2 + val.cosThetaOver2 * val.sinTheta;

                // The vector field derivatives w.r.t. to the local polar coordinates. 
                // The vector components refer to the local cartesian system though.
                var polarGradient = Matrix2by2.CreateFromArray(
                    new double[,] { { a * b_r * c1, a * b * c1_theta }, { a * b_r * c2, a * b * c2_theta } });
                // The vector field derivatives w.r.t. to the local cartesian coordinates.
                displacementGradientMode2 =
                    polarJacobians.TransformVectorFieldDerivativesLocalPolarToLocalCartesian(polarGradient);
            }
        }

        private static Tensor2D ComputeStrainTensor(Matrix2by2 displacementGradient)
        {
            double exx = displacementGradient[0, 0];
            double eyy = displacementGradient[1, 1];
            double exy = 0.5 * (displacementGradient[0, 1] + displacementGradient[1, 0]);
            return new Tensor2D(exx, eyy, exy);
        }

        private static void ComputeStressTensors(CommonValues val, 
            out Tensor2D stressTensorMode1, out Tensor2D stressTensorMode2)
        {
            double coeff = 1.0 / (Math.Sqrt(2.0 * Math.PI) * val.sqrtR);

            double sxxMode1 = coeff * val.cosThetaOver2 * (1.0 - val.sinThetaOver2 * val.sin3ThetaOver2);
            double syyMode1 = coeff * val.cosThetaOver2 * (1.0 + val.sinThetaOver2 * val.sin3ThetaOver2);
            double sxyMode1 = coeff * val.sinThetaOver2 * val.cosThetaOver2 * val.cos3ThetaOver2;

            double sxxMode2 = -coeff * val.sinThetaOver2 * (2.0 + val.cosThetaOver2 * val.cos3ThetaOver2);
            double syyMode2 = sxyMode1;
            double sxyMode2 = sxxMode1;

            stressTensorMode1 = new Tensor2D(sxxMode1, syyMode1, sxyMode1);
            stressTensorMode2 = new Tensor2D(sxxMode2, syyMode2, sxyMode2);
        }

        /// <summary>
        /// A DTO to pass the common values to the various methods, instead of recalculating them. 
        /// Also handles their calculation.
        /// </summary>
        private class CommonValues
        {
            public readonly double sqrtR;
            public readonly double cosTheta;
            public readonly double sinTheta;
            public readonly double cosThetaOver2;
            public readonly double sinThetaOver2;
            public readonly double cos3ThetaOver2;
            public readonly double sin3ThetaOver2;

            public CommonValues(double r, double theta)
            {
                sqrtR = Math.Sqrt(r);
                cosTheta = Math.Cos(theta);
                sinTheta = Math.Sin(theta);
                cosThetaOver2 = Math.Cos(theta / 2.0);
                sinThetaOver2 = Math.Sin(theta / 2.0);
                cos3ThetaOver2 = Math.Cos(3.0 * theta / 2.0);
                sin3ThetaOver2 = Math.Sin(3.0 * theta / 2.0);
            }
        }
    }
}
