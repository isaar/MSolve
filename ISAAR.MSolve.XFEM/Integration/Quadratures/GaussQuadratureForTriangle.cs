using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Integration.Points;

namespace ISAAR.MSolve.XFEM.Integration.Quadratures
{
    sealed class GaussQuadratureForTriangle: IStandardQuadrature2D
    {
        public static readonly GaussQuadratureForTriangle Order1Point1 = new GaussQuadratureForTriangle(
            new GaussPoint2D(1.0/3,     1.0/3,      1.0/2));

        public static readonly GaussQuadratureForTriangle Order2Points3 = new GaussQuadratureForTriangle(
            new GaussPoint2D(1.0 / 6,   1.0 / 6,    1.0 / 6),
            new GaussPoint2D(2.0 / 3,   1.0 / 6,    1.0 / 6),
            new GaussPoint2D(1.0 / 6,   2.0 / 3,    1.0 / 6));

        #region  Alternative quadrature(Order=2, Points=3)
        public static readonly GaussQuadratureForTriangle Order2Points3sAlternative = new GaussQuadratureForTriangle(
            new GaussPoint2D(0.0,       1.0 / 2,    1.0 / 6),
            new GaussPoint2D(1.0 / 2,   0.0,        1.0 / 6),
            new GaussPoint2D(1.0 / 2,   1.0 / 2,    1.0 / 6));
        #endregion

        public static readonly GaussQuadratureForTriangle Order3Points4 = new GaussQuadratureForTriangle(
            new GaussPoint2D(1.0 / 3,   1.0 / 3,    -27.0 / 96),
            new GaussPoint2D(1.0 / 5,   1.0 / 5,    25.0 / 96),
            new GaussPoint2D(1.0 / 5,   3.0 / 5,    25.0 / 96),
            new GaussPoint2D(3.0 / 5,   1.0 / 5,    25.0 / 96));

        public static readonly GaussQuadratureForTriangle Order4Points7 = new GaussQuadratureForTriangle(
            new GaussPoint2D(0.1012865073235, 0.1012865073235, 0.1259391805448 / 2),
            new GaussPoint2D(0.7974269853531, 0.1012865073235, 0.1259391805448 / 2),
            new GaussPoint2D(0.1012865073235, 0.7974269853531, 0.1259391805448 / 2),
            new GaussPoint2D(0.4701420641051, 0.0597158717898, 0.1323941527885 / 2),
            new GaussPoint2D(0.4701420641051, 0.4701420641051, 0.1323941527885 / 2),
            new GaussPoint2D(0.0597158717898, 0.4701420641051, 0.1323941527885 / 2),
            new GaussPoint2D(0.3333333333333, 0.3333333333333, 0.1125));

        #region Alternative quadrature (Order=4, Points=7), from Solverize. Not sure if it is correct.
        public static readonly GaussQuadratureForTriangle Order4Points7Alternative = new GaussQuadratureForTriangle(
            new GaussPoint2D(0.0, 0.0, 1.0 / 40),
            new GaussPoint2D(1.0 / 2, 0.0, 1.0 / 15),
            new GaussPoint2D(1.0, 0.0, 1.0 / 40),
            new GaussPoint2D(1.0 / 2, 1.0 / 2, 1.0 / 15),
            new GaussPoint2D(0.0, 1.0, 1.0 / 40),
            new GaussPoint2D(0.0, 1.0 / 2, 1.0 / 15),
            new GaussPoint2D(1.0 / 3, 1.0 / 3, 9.0 / 40));
        #endregion

        public static readonly GaussQuadratureForTriangle Order7Points13 = new GaussQuadratureForTriangle(
            new GaussPoint2D(0.0651301029022, 0.0651301029022, 0.0533472356088 / 2.0),
            new GaussPoint2D(0.8697397941956, 0.0651301029022, 0.0533472356088 / 2.0),
            new GaussPoint2D(0.0651301029022, 0.8697397941956, 0.0533472356088 / 2.0),
            new GaussPoint2D(0.3128654960049, 0.0486903154253, 0.0771137608903 / 2.0),
            new GaussPoint2D(0.6384441885698, 0.3128654960049, 0.0771137608903 / 2.0),
            new GaussPoint2D(0.0486903154253, 0.6384441885698, 0.0771137608903 / 2.0),
            new GaussPoint2D(0.6384441885698, 0.0486903154253, 0.0771137608903 / 2.0),
            new GaussPoint2D(0.3128654960049, 0.6384441885698, 0.0771137608903 / 2.0),
            new GaussPoint2D(0.0486903154253, 0.3128654960049, 0.0771137608903 / 2.0),
            new GaussPoint2D(0.2603459660790, 0.2603459660790, 0.1756152576332 / 2.0),
            new GaussPoint2D(0.4793080678419, 0.2603459660790, 0.1756152576332 / 2.0),
            new GaussPoint2D(0.2603459660790, 0.4793080678419, 0.1756152576332 / 2.0),
            new GaussPoint2D(0.3333333333333, 0.3333333333333, -0.1495700444677 / 2.0));


        public IReadOnlyList<GaussPoint2D> IntegrationPoints { get; }

        private GaussQuadratureForTriangle(params GaussPoint2D[] points)
        {
            this.IntegrationPoints = new List<GaussPoint2D>(points);
        }
    }
}
