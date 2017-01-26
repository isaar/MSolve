using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Integration.GaussPoints
{
    sealed class IntegrationRule1D
    {
        public static readonly IntegrationRule1D Order1 = new IntegrationRule1D(
            new GaussPoint1D(0.0, 2.0));

        public static readonly IntegrationRule1D Order2 = new IntegrationRule1D(
            new GaussPoint1D(-0.5773502691896, 1.0), 
            new GaussPoint1D(0.5773502691896, 1.0));

        public static readonly IntegrationRule1D Order3 = new IntegrationRule1D(
            new GaussPoint1D(-0.7745966692415, 0.5555555555556),
            new GaussPoint1D(0, 0.8888888888889),
            new GaussPoint1D(0.7745966692415, 0.5555555555556));

        public static readonly IntegrationRule1D Order4 = new IntegrationRule1D(
            new GaussPoint1D(-0.86113631159416, 0.3478548451375),
            new GaussPoint1D(-0.3399810435849, 0.6521451548625),
            new GaussPoint1D(0.3399810435849, 0.6521451548625),
            new GaussPoint1D(0.86113631159416, 0.3478548451375));


        public IReadOnlyList<GaussPoint1D> Points { get; }

        private IntegrationRule1D(params GaussPoint1D[] points)
        {
            this.Points = new List<GaussPoint1D>(points);
        }

    }
}
