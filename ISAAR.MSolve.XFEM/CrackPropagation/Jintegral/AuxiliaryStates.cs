using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Items.CrackTip;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.LinearAlgebra;
using ISAAR.MSolve.XFEM.Materials.Crack;
using ISAAR.MSolve.XFEM.Tensors;

namespace ISAAR.MSolve.XFEM.CrackPropagation.Jintegral
{
    /// <summary>
    /// All quantities are with respect to the local cartesian system of the crack tip.
    /// </summary>
    class AuxiliaryStates
    {
        /// <summary>
        /// The derivatives of the displacement field of an imaginary pure Mode I (opening crack extension) state 
        /// represented in the local cartesian system of the crack tip. The differentation is also w.r.t. the tip local
        /// cartesian coordinates. Matrix dimensions: 2x2.
        /// </summary>
        public DenseMatrix DisplacementGradientMode1 { get; }

        /// <summary>
        /// The derivatives of the displacement field of an imaginary pure Mode II (sliding crack extension) state 
        /// represented in the local cartesian system of the crack tip. The differentation is also w.r.t. the tip local 
        /// cartesian coordinates. Matrix dimensions: 2x2.
        /// </summary>
        public DenseMatrix DisplacementGradientMode2 { get; }

        /// <summary>
        /// The strain tensor of an imaginary pure Mode I (opening crack extension) state represented in the local
        /// cartesian system of the crack tip. Tensor dimensions: 2x2.
        /// </summary>
        public Tensor2D StrainTensorMode1 { get; }

        /// <summary>
        /// The strain tensor of an imaginary pure Mode IΙ (sliding crack extension) state represented in the local
        /// cartesian system of the crack tip. Tensor dimensions: 2x2.
        /// </summary>
        public Tensor2D StrainTensorMode2 { get; }

        /// <summary>
        /// The stress tensor of an imaginary pure Mode I (opening crack extension) state represented in the local
        /// cartesian system of the crack tip. Tensor dimensions: 2x2.
        /// </summary>
        public Tensor2D StressTensorMode1 { get; }

        /// <summary>
        /// The stress tensor of an imaginary pure Mode IΙ (sliding crack extension) state represented in the local
        /// cartesian system of the crack tip. Tensor dimensions: 2x2.
        /// </summary>
        public Tensor2D StressTensorMode2 { get; }

        public static AuxiliaryStates ComputeLocalCartesianTensorsOfAuxiliaryStates(CrackMaterial2D materialAtTip, 
            TipCoordinateSystem tipCoordinateSystem, ICartesianPoint2D localCartesianCoordinatesOfGaussPoint)
        {
            return new AuxiliaryStates(materialAtTip, tipCoordinateSystem, localCartesianCoordinatesOfGaussPoint);
        }

        private AuxiliaryStates(CrackMaterial2D materialAtTip, TipCoordinateSystem tipCoordinateSystem, 
            ICartesianPoint2D localCartesianCoordinatesOfGaussPoint)
        {
            // Common calculations
            PolarPoint2D point = 
                tipCoordinateSystem.TransformPointLocalCartesianToLocalPolar(localCartesianCoordinatesOfGaussPoint);
            var commonValues = new CommonValues(point.R, point.Theta);

            // Displacement field derivatives
            TipJacobians polarJacobians = tipCoordinateSystem.CalculateJacobiansAt(point);
            Tuple<DenseMatrix, DenseMatrix> gradients = 
                ComputeDisplacementDerivatives(polarJacobians, materialAtTip, commonValues);
            DisplacementGradientMode1 = gradients.Item1;
            DisplacementGradientMode2 = gradients.Item2;

            // Strains
            StrainTensorMode1 = ComputeStrainTensor(DisplacementGradientMode1);
            StrainTensorMode2 = ComputeStrainTensor(DisplacementGradientMode2);

            // Stresses
            Tuple<Tensor2D, Tensor2D> stresses = ComputeStressTensors(commonValues);
            StressTensorMode1 = stresses.Item1;
            StressTensorMode2 = stresses.Item2;
        }

        private static Tuple<DenseMatrix, DenseMatrix> ComputeDisplacementDerivatives(TipJacobians polarJacobians, 
            CrackMaterial2D materialAtTip, CommonValues val)
        {
            DenseMatrix displacementGradientMode1, displacementGradientMode2;

            // Temporary values and derivatives of the differentiated quantities. See documentation for their derivation.
            double k = materialAtTip.KolosovCoefficient;
            double a = (1.0 + materialAtTip.PoissonRatio) / (materialAtTip.YoungModulus * Math.Sqrt(2.0 * Math.PI));
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
                DenseMatrix polarGradient = new DenseMatrix(new double[,] { { a * b_r * c1, a * b * c1_theta },
                                                                        { a * b_r * c2, a * b * c2_theta } });
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
                DenseMatrix polarGradient = new DenseMatrix(new double[,] { { a * b_r * c1, a * b * c1_theta },
                                                                        { a * b_r * c2, a * b * c2_theta } });
                // The vector field derivatives w.r.t. to the local cartesian coordinates.
                displacementGradientMode2 =
                    polarJacobians.TransformVectorFieldDerivativesLocalPolarToLocalCartesian(polarGradient);
            }

            return new Tuple<DenseMatrix, DenseMatrix>(displacementGradientMode1, displacementGradientMode2);
        }

        private static Tensor2D ComputeStrainTensor(DenseMatrix displacementGradient)
        {
            double exx = displacementGradient[0, 0];
            double eyy = displacementGradient[1, 1];
            double exy = 0.5 * (displacementGradient[0, 1] + displacementGradient[1, 0]);
            return new Tensor2D(exx, eyy, exy);
        }

        private static Tuple<Tensor2D, Tensor2D> ComputeStressTensors(CommonValues val)
        {
            double coeff = 1.0 / (Math.Sqrt(2.0 * Math.PI) * val.sqrtR);

            double sxxMode1 = coeff * val.cosThetaOver2 * (1.0 - val.sinThetaOver2 * val.sin3ThetaOver2);
            double syyMode1 = coeff * val.cosThetaOver2 * (1.0 + val.sinThetaOver2 * val.sin3ThetaOver2);
            double sxyMode1 = coeff * val.sinThetaOver2 * val.cosThetaOver2 * val.cos3ThetaOver2;

            double sxxMode2 = -coeff * val.sinThetaOver2 * (2.0 + val.cosThetaOver2 * val.cos3ThetaOver2);
            double syyMode2 = sxyMode1;
            double sxyMode2 = sxxMode1;

            return new Tuple<Tensor2D, Tensor2D>(
                new Tensor2D(sxxMode1, syyMode1, sxyMode1), new Tensor2D(sxxMode2, syyMode2, sxyMode2));
        }

        /// <summary>
        /// A DTO to pass the common values to the various methods, instead of recalculating them
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
