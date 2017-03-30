using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Items.CrackTip;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Enrichments.Functions
{
    // TODO: The cartesian global to cartesian local to polar system transformation and the jacobian should be 
    // calculated only once per gauss point for all functions of the same enrichment item.
    static class IsotropicBrittleTipFunctions2DAlternative
    {
        public class Func1 : IEnrichmentFunction2D
        {
            private readonly CrackTip2D enrichmentItem;

            public Func1(CrackTip2D enrichmentItem)
            {
                this.enrichmentItem = enrichmentItem;
            }

            public double EvalueAt(ICartesianPoint2D cartesianPoint)
            {
                //This section needs to be calculated only once per gauss point for all enrichment items of the same item
                PolarPoint2D polarCoordinates = 
                    enrichmentItem.TipSystem.TransformPointGlobalCartesianToLocalPolar(cartesianPoint);

                return Math.Sqrt(polarCoordinates.R) * Math.Sin(polarCoordinates.Theta / 2.0);
            }

            public EvaluatedFunction2D EvaluateAllAt(ICartesianPoint2D cartesianPoint)
            {
                //This section needs to be calculated only once per gauss point for all enrichment items of the same item
                PolarPoint2D polarCoordinates =
                    enrichmentItem.TipSystem.TransformPointGlobalCartesianToLocalPolar(cartesianPoint);
                TipJacobians jacobian = enrichmentItem.TipSystem.CalculateJacobiansAt(polarCoordinates);

                double sqrtR = Math.Sqrt(polarCoordinates.R);
                double cosThetaHalf = Math.Cos(polarCoordinates.Theta / 2.0);
                double sinThetaHalf = Math.Sin(polarCoordinates.Theta / 2.0);

                double value = sqrtR * sinThetaHalf;
                double derivativeR = 0.5 / sqrtR * sinThetaHalf;
                double derivativeTheta = 0.5 * sqrtR * cosThetaHalf;

                return new EvaluatedFunction2D(value, 
                    jacobian.TransformScalarFieldDerivativesLocalPolarToGlobalCartesian(derivativeR, derivativeTheta));
            }
        }

        public class Func2 : IEnrichmentFunction2D
        {
            private readonly CrackTip2D enrichmentItem;

            public Func2(CrackTip2D enrichmentItem)
            {
                this.enrichmentItem = enrichmentItem;
            }

            public double EvalueAt(ICartesianPoint2D cartesianPoint)
            {
                //This section needs to be calculated only once per gauss point for all enrichment items of the same item
                PolarPoint2D polarCoordinates =
                    enrichmentItem.TipSystem.TransformPointGlobalCartesianToLocalPolar(cartesianPoint);

                return Math.Sqrt(polarCoordinates.R) * Math.Cos(polarCoordinates.Theta / 2.0);
            }

            public EvaluatedFunction2D EvaluateAllAt(ICartesianPoint2D cartesianPoint)
            {
                //This section needs to be calculated only once per gauss point for all enrichment items of the same item
                PolarPoint2D polarCoordinates =
                    enrichmentItem.TipSystem.TransformPointGlobalCartesianToLocalPolar(cartesianPoint);
                TipJacobians jacobian = enrichmentItem.TipSystem.CalculateJacobiansAt(polarCoordinates);

                double sqrtR = Math.Sqrt(polarCoordinates.R);
                double cosThetaHalf = Math.Cos(polarCoordinates.Theta / 2.0);
                double sinThetaHalf = Math.Sin(polarCoordinates.Theta / 2.0);
                

                double value = sqrtR * cosThetaHalf;
                double derivativeR = 0.5 / sqrtR * cosThetaHalf;
                double derivativeTheta = -0.5 * sqrtR * sinThetaHalf;

                return new EvaluatedFunction2D(value,
                    jacobian.TransformScalarFieldDerivativesLocalPolarToGlobalCartesian(derivativeR, derivativeTheta));
            }
        }

        public class Func3 : IEnrichmentFunction2D
        {
            private readonly CrackTip2D enrichmentItem;

            public Func3(CrackTip2D enrichmentItem)
            {
                this.enrichmentItem = enrichmentItem;
            }

            public double EvalueAt(ICartesianPoint2D cartesianPoint)
            {
                //This section needs to be calculated only once per gauss point for all enrichment items of the same item
                PolarPoint2D polarCoordinates =
                    enrichmentItem.TipSystem.TransformPointGlobalCartesianToLocalPolar(cartesianPoint);

                return Math.Sqrt(polarCoordinates.R) * Math.Sin(polarCoordinates.Theta / 2.0) 
                    * Math.Cos(polarCoordinates.Theta);
            }

            public EvaluatedFunction2D EvaluateAllAt(ICartesianPoint2D cartesianPoint)
            {
                //This section needs to be calculated only once per gauss point for all enrichment items of the same item
                PolarPoint2D polarCoordinates =
                    enrichmentItem.TipSystem.TransformPointGlobalCartesianToLocalPolar(cartesianPoint);
                TipJacobians jacobian = enrichmentItem.TipSystem.CalculateJacobiansAt(polarCoordinates);

                double sqrtR = Math.Sqrt(polarCoordinates.R);
                double cosTheta = Math.Cos(polarCoordinates.Theta);
                double sinTheta = Math.Sin(polarCoordinates.Theta);
                double cosThetaHalf = Math.Cos(polarCoordinates.Theta / 2.0);
                double sinThetaHalf = Math.Sin(polarCoordinates.Theta / 2.0);

                double value = sqrtR * sinThetaHalf * cosTheta;
                double derivativeR = 0.5 / sqrtR * sinThetaHalf * cosTheta;
                double derivativeTheta = sqrtR * (0.5 * cosThetaHalf * cosTheta - sinThetaHalf * sinTheta);

                return new EvaluatedFunction2D(value,
                    jacobian.TransformScalarFieldDerivativesLocalPolarToGlobalCartesian(derivativeR, derivativeTheta));
            }
        }

        public class Func4 : IEnrichmentFunction2D
        {
            private readonly CrackTip2D enrichmentItem;

            public Func4(CrackTip2D enrichmentItem)
            {
                this.enrichmentItem = enrichmentItem;
            }

            public double EvalueAt(ICartesianPoint2D cartesianPoint)
            {
                //This section needs to be calculated only once per gauss point for all enrichment items of the same item
                PolarPoint2D polarCoordinates =
                    enrichmentItem.TipSystem.TransformPointGlobalCartesianToLocalPolar(cartesianPoint);

                return Math.Sqrt(polarCoordinates.R) * Math.Cos(polarCoordinates.Theta / 2.0)
                    * Math.Cos(polarCoordinates.Theta);
            }

            public EvaluatedFunction2D EvaluateAllAt(ICartesianPoint2D cartesianPoint)
            {
                //This section needs to be calculated only once per gauss point for all enrichment items of the same item
                PolarPoint2D polarCoordinates =
                    enrichmentItem.TipSystem.TransformPointGlobalCartesianToLocalPolar(cartesianPoint);
                TipJacobians jacobian = enrichmentItem.TipSystem.CalculateJacobiansAt(polarCoordinates);

                double sqrtR = Math.Sqrt(polarCoordinates.R);
                double cosTheta = Math.Cos(polarCoordinates.Theta);
                double sinTheta = Math.Sin(polarCoordinates.Theta);
                double cosThetaHalf = Math.Cos(polarCoordinates.Theta / 2.0);
                double sinThetaHalf = Math.Sin(polarCoordinates.Theta / 2.0);

                double value = sqrtR * cosThetaHalf * cosTheta;
                double derivativeR = 0.5 / sqrtR * cosThetaHalf * cosTheta;
                double derivativeTheta = sqrtR * (-0.5 * sinThetaHalf * cosTheta - cosThetaHalf * sinTheta);

                return new EvaluatedFunction2D(value,
                    jacobian.TransformScalarFieldDerivativesLocalPolarToGlobalCartesian(derivativeR, derivativeTheta));
            }
        }
    }
}
