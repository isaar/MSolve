using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Utilities;

//TODO: Replace double[] with Vector
namespace ISAAR.MSolve.XFEM.Enrichments.Functions
{
    // TODO: The cartesian global to cartesian local to polar system transformation and the jacobian should be 
    // calculated only once per gauss point for all functions of the same enrichment item.
    static class IsotropicBrittleTipFunctions2D
    {
        public class Func1 : ITipFunction
        {
            public double EvaluateAt(PolarPoint2D point)
            {
                return Math.Sqrt(point.R) * Math.Sin(point.Theta / 2.0);
            }

            public EvaluatedFunction2D EvaluateAllAt(PolarPoint2D point, TipJacobians jacobian)
            {
                double sqrtR = Math.Sqrt(point.R);
                double cosThetaHalf = Math.Cos(point.Theta / 2.0);
                double sinThetaHalf = Math.Sin(point.Theta / 2.0);

                double value = sqrtR * sinThetaHalf;
                var gradientPolar = Vector2.Create(0.5 / sqrtR * sinThetaHalf, 0.5 * sqrtR * cosThetaHalf);

                Vector2 gradientGlobal = jacobian.TransformScalarFieldDerivativesLocalPolarToGlobalCartesian(gradientPolar);
                return new EvaluatedFunction2D(value, gradientGlobal);
            }

            public override string ToString()
            {
                return "TipIsotropicBrittle1";
            }
        }

        public class Func2 : ITipFunction
        {
            public double EvaluateAt(PolarPoint2D point)
            {
                return Math.Sqrt(point.R) * Math.Cos(point.Theta / 2.0);
            }

            public EvaluatedFunction2D EvaluateAllAt(PolarPoint2D point, TipJacobians jacobian)
            {
                double sqrtR = Math.Sqrt(point.R);
                double cosThetaHalf = Math.Cos(point.Theta / 2.0);
                double sinThetaHalf = Math.Sin(point.Theta / 2.0);

                double value = sqrtR * cosThetaHalf;
                var gradientPolar = Vector2.Create(0.5 / sqrtR * cosThetaHalf, -0.5 * sqrtR * sinThetaHalf);

                Vector2 gradientGlobal = jacobian.TransformScalarFieldDerivativesLocalPolarToGlobalCartesian(gradientPolar);
                return new EvaluatedFunction2D(value, gradientGlobal);
            }

            public override string ToString()
            {
                return "TipIsotropicBrittle2";
            }
        }

        public class Func3 : ITipFunction
        {
            public double EvaluateAt(PolarPoint2D point)
            {
                return Math.Sqrt(point.R) * Math.Sin(point.Theta / 2.0) * Math.Sin(point.Theta);
            }

            public EvaluatedFunction2D EvaluateAllAt(PolarPoint2D point, TipJacobians jacobian)
            {
                double sqrtR = Math.Sqrt(point.R);
                double cosTheta = Math.Cos(point.Theta);
                double sinTheta = Math.Sin(point.Theta);
                double cosThetaHalf = Math.Cos(point.Theta / 2.0);
                double sinThetaHalf = Math.Sin(point.Theta / 2.0);

                double value = sqrtR * sinThetaHalf * sinTheta;
                var gradientPolar = Vector2.Create(0.5 / sqrtR * sinThetaHalf * sinTheta,
                    sqrtR * (0.5 * cosThetaHalf * sinTheta + sinThetaHalf * cosTheta) );

                Vector2 gradientGlobal = jacobian.TransformScalarFieldDerivativesLocalPolarToGlobalCartesian(gradientPolar);
                return new EvaluatedFunction2D(value, gradientGlobal);
            }

            public override string ToString()
            {
                return "TipIsotropicBrittle3";
            }
        }

        public class Func4 : ITipFunction
        {
            public double EvaluateAt(PolarPoint2D point)
            {
                return Math.Sqrt(point.R) * Math.Cos(point.Theta / 2.0) * Math.Sin(point.Theta);
            }

            public EvaluatedFunction2D EvaluateAllAt(PolarPoint2D point, TipJacobians jacobian)
            {
                double sqrtR = Math.Sqrt(point.R);
                double cosTheta = Math.Cos(point.Theta);
                double sinTheta = Math.Sin(point.Theta);
                double cosThetaHalf = Math.Cos(point.Theta / 2.0);
                double sinThetaHalf = Math.Sin(point.Theta / 2.0);

                double value = sqrtR * cosThetaHalf * sinTheta;
                var gradientPolar = Vector2.Create(0.5 / sqrtR * cosThetaHalf * sinTheta,
                    sqrtR * (-0.5 * sinThetaHalf * sinTheta + cosThetaHalf * cosTheta) );

                Vector2 gradientGlobal = jacobian.TransformScalarFieldDerivativesLocalPolarToGlobalCartesian(gradientPolar);
                return new EvaluatedFunction2D(value, gradientGlobal);
            }

            public override string ToString()
            {
                return "TipIsotropicBrittle4";
            }
        }
    }
}