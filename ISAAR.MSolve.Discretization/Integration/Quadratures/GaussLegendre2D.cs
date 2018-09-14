using ISAAR.MSolve.Discretization.Integration.Points;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Discretization.Integration.Quadratures
{
    /// <summary>
    /// Enum class with the 2D Gauss-Legendre integration rules of varying orders. Tensor product of 
    /// <see cref="GaussLegendre1D"/> rules along each axis.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public sealed class GaussLegendre2D : IQuadrature2D
    {
        public static readonly GaussLegendre2D Order1x1 = new GaussLegendre2D(GaussLegendre1D.Order1, GaussLegendre1D.Order1);
        public static readonly GaussLegendre2D Order2x2 = new GaussLegendre2D(GaussLegendre1D.Order2, GaussLegendre1D.Order2);
        public static readonly GaussLegendre2D Order3x3 = new GaussLegendre2D(GaussLegendre1D.Order3, GaussLegendre1D.Order3);
        public static readonly GaussLegendre2D Order4x4 = new GaussLegendre2D(GaussLegendre1D.Order4, GaussLegendre1D.Order4);
        public static readonly GaussLegendre2D Order5x5 = new GaussLegendre2D(GaussLegendre1D.Order5, GaussLegendre1D.Order5);
        public static readonly GaussLegendre2D Order6x6 = new GaussLegendre2D(GaussLegendre1D.Order6, GaussLegendre1D.Order6);
        public static readonly GaussLegendre2D Order7x7 = new GaussLegendre2D(GaussLegendre1D.Order7, GaussLegendre1D.Order7);
        public static readonly GaussLegendre2D Order8x8 = new GaussLegendre2D(GaussLegendre1D.Order8, GaussLegendre1D.Order8);
        public static readonly GaussLegendre2D Order9x9 = new GaussLegendre2D(GaussLegendre1D.Order9, GaussLegendre1D.Order9);
        public static readonly GaussLegendre2D Order10x10 = new GaussLegendre2D(GaussLegendre1D.Order10, GaussLegendre1D.Order10);

        private GaussLegendre2D(GaussLegendre1D ruleXi, GaussLegendre1D ruleEta)
        {
            // Combine the integration rules of each axis: 
            // WARNING: Do not change their order (Xi major, Eta minor). Other classes, such as ExtrapolationGaussLegendre2x2 
            //          depend on it.
            var points2D = new List<GaussPoint2D>();
            foreach (var pointEta in ruleEta.IntegrationPoints)
            {
                foreach (var pointXi in ruleXi.IntegrationPoints)
                {
                    points2D.Add(new GaussPoint2D(pointXi.Xi, pointEta.Xi, pointXi.Weight * pointEta.Weight));
                }
            }
            this.IntegrationPoints = points2D;
        }

        /// <summary>
        /// The integration points are sorted based on an order strictly defined for each quadrature.
        /// </summary>
        public IReadOnlyList<GaussPoint2D> IntegrationPoints { get; }
    }
}
