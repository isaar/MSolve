using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Numerical.Commons;

//TODO: A thread safe Table3D is needed with an atomic method: 
//      GetOrAddNew(int orderXi, int orderEta, int orderZeta Func<int, int, int, GaussLegendre3D> createFunc). This quadrature  
//      class should also be thread safe. What about distributed systems?
namespace ISAAR.MSolve.Discretization.Integration.Quadratures
{
    /// <summary>
    /// Contains the 3D Gauss-Legendre integration rules of varying orders. These are created as tensor products of 
    /// <see cref="GaussLegendre1D"/> rules along each axis. Each quadrature is a singleton that can be accessed through this 
    /// class.
    /// Authors: Dimitris Tsapetis, Serafeim Bakalakos
    /// </summary>
    public sealed class GaussLegendre3D : IQuadrature3D
    {
        private const int initialCapacity = 10;

        private static readonly Table3D<int, int, int, GaussLegendre3D> quadratures =
            new Table3D<int, int, int, GaussLegendre3D>(initialCapacity);

        public static GaussLegendre3D GetQuadratureWithOrder(int orderXi, int orderEta, int orderZeta)
        {
            //TODO: these should be thread safe and atomic.
            bool exists = quadratures.TryGetValue(orderXi, orderEta, orderZeta, out GaussLegendre3D quadrature);
            if (!exists)
            {
                quadrature = new GaussLegendre3D(orderXi, orderEta, orderZeta);
                quadratures[orderXi, orderEta, orderZeta] = quadrature;
            }
            return quadrature;
        }

        private GaussLegendre3D(int orderXi, int orderEta, int orderZeta)
        {
            GaussLegendre1D quadratureXi = GaussLegendre1D.GetQuadratureWithOrder(orderXi);
            GaussLegendre1D quadratureEta = GaussLegendre1D.GetQuadratureWithOrder(orderEta);
            GaussLegendre1D quadratureZeta = GaussLegendre1D.GetQuadratureWithOrder(orderZeta);
            var points3D = new List<GaussPoint3D>();

            // Combine the integration rules of each axis. The order is Xi minor - Eta middle - Zeta major
            // WARNING: Do not change their order. Other classes, such as the ones that implement extrapolations, depend on it.
            foreach (var pointZeta in quadratureZeta.IntegrationPoints)
            {
                foreach (var pointEta in quadratureEta.IntegrationPoints)
                {
                    foreach (var pointXi in quadratureXi.IntegrationPoints)
                    {
                        points3D.Add(new GaussPoint3D(pointXi.Xi, pointEta.Xi, pointZeta.Xi,
                            pointXi.Weight * pointEta.Weight * pointZeta.Weight));
                    }
                }
            }

            this.IntegrationPoints = points3D;
        }

        /// <summary>
        /// The integration points are sorted based on an order strictly defined for each quadrature.
        /// </summary>
        public IReadOnlyList<GaussPoint3D> IntegrationPoints { get; }
    }
}