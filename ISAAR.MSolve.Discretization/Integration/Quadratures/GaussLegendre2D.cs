using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Numerical.Commons;

//TODO: A thread safe Table is needed with an atomic method: 
//      GetOrAddNew(int orderXi, int orderEta, Func<int, int, GaussLegendre2D> createFunc). This quadrature class should also 
//      be thread safe. What about distributed systems?
namespace ISAAR.MSolve.Discretization.Integration.Quadratures
{
    /// <summary>
    /// Contains the 2D Gauss-Legendre integration rules of varying orders. These are created as tensor products of 
    /// <see cref="GaussLegendre1D"/> rules along each axis. Each quadrature is a singleton that can be accessed through this 
    /// class.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public sealed class GaussLegendre2D : IQuadrature2D
    {
        private const int initialCapacity = 10;

        private static readonly Table<int, int, GaussLegendre2D> quadratures = 
            new Table<int, int, GaussLegendre2D>(initialCapacity); 

        public static GaussLegendre2D GetQuadratureWithOrder(int orderXi, int orderEta)
        {
            //TODO: these should be thread safe and atomic.
            bool exists = quadratures.TryGetValue(orderXi, orderEta, out GaussLegendre2D quadrature);
            if (!exists)
            {
                quadrature = new GaussLegendre2D(orderXi, orderEta);
                quadratures[orderXi, orderEta] = quadrature;
            }
            return quadrature;
        }

        private GaussLegendre2D(int orderXi, int orderEta)
        {
            GaussLegendre1D quadratureXi = GaussLegendre1D.GetQuadratureWithOrder(orderXi);
            GaussLegendre1D quadratureEta = GaussLegendre1D.GetQuadratureWithOrder(orderEta);
            var points2D = new List<GaussPoint2D>();

            // Combine the integration rules of each axis. The order is Xi minor - Eta major 
            // WARNING: Do not change their order. Other classes, such as ExtrapolationGaussLegendre2x2 depend on it.
            foreach (var pointEta in quadratureXi.IntegrationPoints)
            {
                foreach (var pointXi in quadratureEta.IntegrationPoints)
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
