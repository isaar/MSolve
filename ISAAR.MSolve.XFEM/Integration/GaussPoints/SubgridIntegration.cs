using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Integration.GaussPoints
{
    class SubgridIntegration
    {
        private readonly int subgridsPerAxis;
        private readonly IntegrationRule2D gaussQuadrature;

        public SubgridIntegration(int subgridsPerAxis)
        {
            this.subgridsPerAxis = subgridsPerAxis;
            this.gaussQuadrature = IntegrationRule2D.Order2x2;
        }

        public SubgridIntegration(int subgridsPerAxis, IntegrationRule2D gaussQuadrature)
        {
            this.subgridsPerAxis = subgridsPerAxis;
            this.gaussQuadrature = gaussQuadrature;
        }

        public IReadOnlyList<GaussPoint2D> GeneratePoints()
        {
            var points = new List<GaussPoint2D>();
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

                    foreach(var subgridPoint in gaussQuadrature.Points)
                    {
                        // Transformation from the system of the subrectangle to the natural system of the element
                        double naturalXi = subgridPoint.X * (xiMax - xiMin) / 2.0 + (xiMin + xiMax) / 2.0;
                        double naturalEta = subgridPoint.Y * (etaMax - etaMin) / 2.0 + (etaMin + etaMax) / 2.0;
                        double naturalWeight = subgridPoint.Weight * (xiMax - xiMin) / 2.0 * (etaMax - etaMin) / 2.0;
                        points.Add(new GaussPoint2D(naturalXi, naturalEta, naturalWeight));
                    }
                }
            }
            return points;
        }
    }
}
