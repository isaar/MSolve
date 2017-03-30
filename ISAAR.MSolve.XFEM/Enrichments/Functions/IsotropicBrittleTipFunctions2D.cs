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
    static class IsotropicBrittleTipFunctions2D
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

            public Tuple<double, double> EvaluateDerivativesAt(ICartesianPoint2D cartesianPoint)
            {
                throw new NotImplementedException();
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

            public Tuple<double, double> EvaluateDerivativesAt(ICartesianPoint2D cartesianPoint)
            {
                throw new NotImplementedException();
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
                    * Math.Sin(polarCoordinates.Theta);
            }

            public Tuple<double, double> EvaluateDerivativesAt(ICartesianPoint2D cartesianPoint)
            {
                throw new NotImplementedException();
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

                double value = sqrtR * sinThetaHalf * sinTheta;
                double derivativeR = 0.5 / sqrtR * sinThetaHalf * sinTheta;
                double derivativeTheta = sqrtR * (0.5 * cosThetaHalf * sinTheta + sinThetaHalf * cosTheta);

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
                    * Math.Sin(polarCoordinates.Theta);
            }

            public Tuple<double, double> EvaluateDerivativesAt(ICartesianPoint2D cartesianPoint)
            {
                throw new NotImplementedException();
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

                double value = sqrtR * cosThetaHalf * sinTheta;
                double derivativeR = 0.5 / sqrtR * cosThetaHalf * sinTheta;
                double derivativeTheta = sqrtR * (-0.5 * sinThetaHalf * sinTheta + cosThetaHalf * cosTheta);

                return new EvaluatedFunction2D(value,
                    jacobian.TransformScalarFieldDerivativesLocalPolarToGlobalCartesian(derivativeR, derivativeTheta));
            }
        }
    }
}