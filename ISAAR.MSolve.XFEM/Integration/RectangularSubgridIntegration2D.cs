using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration;
using ISAAR.MSolve.Discretization.Integration.Quadratures;

namespace ISAAR.MSolve.XFEM.Integration
{
    /// <summary>
    /// TODO: This rule is actually independent from the element and its elements can be cached, albeit not in a 
    /// static manner. Should I put it with the standard quadratures?
    /// TODO: Ensure this is not used for anything other than Quadrilaterals.
    /// </summary>
    public class RectangularSubgridIntegration2D<TElement> : IIntegrationStrategy2D<TElement>
    {
        private readonly int subgridsPerAxis;
        private readonly GaussLegendre2D  gaussQuadrature;

        public RectangularSubgridIntegration2D(int subgridsPerAxis) : 
            this(subgridsPerAxis, GaussLegendre2D.GetQuadratureWithOrder(2,2))
        {
        }

        public RectangularSubgridIntegration2D(int subgridsPerAxis, GaussLegendre2D gaussQuadrature)
        {
            this.subgridsPerAxis = subgridsPerAxis;
            this.gaussQuadrature = gaussQuadrature;
        }

        public IReadOnlyList<GaussPoint> GenerateIntegrationPoints(TElement element)
        {
            var points = new List<GaussPoint>();
            double length = 2.0 / subgridsPerAxis;
            for (int i = 0; i < subgridsPerAxis; ++i)
            {
                for (int j = 0; j < subgridsPerAxis; ++j)
                {
                    // The borders of the subrectangle
                    double xiMin = -1.0 + length * i;
                    double xiMax = -1.0 + length * (i+1);
                    double etaMin = -1.0 + length * j;
                    double etaMax = -1.0 + length * (j + 1);

                    foreach(var subgridPoint in gaussQuadrature.IntegrationPoints)
                    {
                        // Transformation from the system of the subrectangle to the natural system of the element
                        double naturalXi = subgridPoint.Xi * (xiMax - xiMin) / 2.0 + (xiMin + xiMax) / 2.0;
                        double naturalEta = subgridPoint.Eta * (etaMax - etaMin) / 2.0 + (etaMin + etaMax) / 2.0;
                        double naturalWeight = subgridPoint.Weight * (xiMax - xiMin) / 2.0 * (etaMax - etaMin) / 2.0;
                        points.Add(new GaussPoint(naturalXi, naturalEta, naturalWeight));
                    }
                }
            }
            return points;
        }
    }
}
