using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Integration.Points;

namespace ISAAR.MSolve.XFEM.Integration.Quadratures
{
    using static GaussLegendre1D;

    sealed class GaussLegendre2D : IStandardQuadrature2D
    {
        public static readonly GaussLegendre2D  Order1x1 = new GaussLegendre2D (Order1, Order1);
        public static readonly GaussLegendre2D  Order2x2 = new GaussLegendre2D (Order2, Order2);
        public static readonly GaussLegendre2D  Order3x3 = new GaussLegendre2D (Order3, Order3);
        public static readonly GaussLegendre2D  Order4x4 = new GaussLegendre2D (Order4, Order4);
        public static readonly GaussLegendre2D  Order5x5 = new GaussLegendre2D (Order5, Order5);
        public static readonly GaussLegendre2D  Order6x6 = new GaussLegendre2D (Order6, Order6);
        public static readonly GaussLegendre2D  Order7x7 = new GaussLegendre2D (Order7, Order7);
        public static readonly GaussLegendre2D  Order8x8 = new GaussLegendre2D (Order8, Order8);
        public static readonly GaussLegendre2D  Order9x9 = new GaussLegendre2D (Order9, Order9);
        public static readonly GaussLegendre2D  Order10x10 = new GaussLegendre2D (Order10, Order10);

        public IReadOnlyList<GaussPoint2D> IntegrationPoints { get; }

        private GaussLegendre2D (GaussLegendre1D ruleXi, GaussLegendre1D ruleEta)
        {
            // Combine the integration rules of each axis: 
            var points2D = new List<GaussPoint2D>();
            foreach (var pointXi in ruleXi.IntegrationPoints)
            {
                foreach (var pointEta in ruleEta.IntegrationPoints)
                {
                    points2D.Add(new GaussPoint2D(pointXi.Xi, pointEta.Xi, pointXi.Weight * pointEta.Weight));
                }
            }
            this.IntegrationPoints = points2D;
        }
    }
}
