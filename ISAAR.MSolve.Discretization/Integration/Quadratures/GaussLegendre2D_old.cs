using ISAAR.MSolve.Discretization.Integration.Points;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Discretization.Integration.Quadratures
{
    /// <summary>
    /// Enum class with the 2D Gauss-Legendre integration rules of varying orders. Tensor product of 
    /// <see cref="GaussLegendre1D_old"/> rules along each axis.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public sealed class GaussLegendre2D_old : IQuadrature2D
    {
        public static readonly GaussLegendre2D_old Order1x1 = new GaussLegendre2D_old(GaussLegendre1D_old.Order1, GaussLegendre1D_old.Order1);
        public static readonly GaussLegendre2D_old Order2x2 = new GaussLegendre2D_old(GaussLegendre1D_old.Order2, GaussLegendre1D_old.Order2);
        public static readonly GaussLegendre2D_old Order3x3 = new GaussLegendre2D_old(GaussLegendre1D_old.Order3, GaussLegendre1D_old.Order3);
        public static readonly GaussLegendre2D_old Order4x4 = new GaussLegendre2D_old(GaussLegendre1D_old.Order4, GaussLegendre1D_old.Order4);
        public static readonly GaussLegendre2D_old Order5x5 = new GaussLegendre2D_old(GaussLegendre1D_old.Order5, GaussLegendre1D_old.Order5);
        public static readonly GaussLegendre2D_old Order6x6 = new GaussLegendre2D_old(GaussLegendre1D_old.Order6, GaussLegendre1D_old.Order6);
        public static readonly GaussLegendre2D_old Order7x7 = new GaussLegendre2D_old(GaussLegendre1D_old.Order7, GaussLegendre1D_old.Order7);
        public static readonly GaussLegendre2D_old Order8x8 = new GaussLegendre2D_old(GaussLegendre1D_old.Order8, GaussLegendre1D_old.Order8);
        public static readonly GaussLegendre2D_old Order9x9 = new GaussLegendre2D_old(GaussLegendre1D_old.Order9, GaussLegendre1D_old.Order9);
        public static readonly GaussLegendre2D_old Order10x10 = new GaussLegendre2D_old(GaussLegendre1D_old.Order10, GaussLegendre1D_old.Order10);

        private GaussLegendre2D_old(GaussLegendre1D_old ruleXi, GaussLegendre1D_old ruleEta)
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
