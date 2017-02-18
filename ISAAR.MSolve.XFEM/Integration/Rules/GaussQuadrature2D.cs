using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Integration.Rules
{
    using Points;
    using static GaussQuadrature1D;

    class GaussQuadrature2D: IIntegrationRule2D
    {
        public static readonly GaussQuadrature2D Order1x1 = new GaussQuadrature2D(Order1, Order1);
        public static readonly GaussQuadrature2D Order2x2 = new GaussQuadrature2D(Order2, Order2);
        public static readonly GaussQuadrature2D Order3x3 = new GaussQuadrature2D(Order3, Order3);
        public static readonly GaussQuadrature2D Order4x4 = new GaussQuadrature2D(Order4, Order4);
        public static readonly GaussQuadrature2D Order5x5 = new GaussQuadrature2D(Order5, Order5);
        public static readonly GaussQuadrature2D Order6x6 = new GaussQuadrature2D(Order6, Order6);
        public static readonly GaussQuadrature2D Order7x7 = new GaussQuadrature2D(Order7, Order7);
        public static readonly GaussQuadrature2D Order8x8 = new GaussQuadrature2D(Order8, Order8);
        public static readonly GaussQuadrature2D Order9x9 = new GaussQuadrature2D(Order9, Order9);
        public static readonly GaussQuadrature2D Order10x10 = new GaussQuadrature2D(Order10, Order10);

        public readonly IReadOnlyList<GaussPoint2D> integrationPoints;

        private GaussQuadrature2D(GaussQuadrature1D ruleXi, GaussQuadrature1D ruleEta)
        {
            // Combine the integration rules of each axis: 
            var points2D = new List<GaussPoint2D>();
            foreach (var pointXi in ruleXi.GenerateIntegrationPoints())
            {
                foreach (var pointEta in ruleEta.GenerateIntegrationPoints())
                {
                    points2D.Add(new GaussPoint2D(pointXi.Xi, pointEta.Xi, pointXi.Weight * pointEta.Weight));
                }
            }
            this.integrationPoints = points2D;
        }

        public IReadOnlyList<GaussPoint2D> GenerateIntegrationPoints() { return integrationPoints; }
    }
}
