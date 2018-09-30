using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Points;

namespace ISAAR.MSolve.Discretization.Integration.Quadratures
{
    /// <summary>
    /// Enum class with the 2D integration rules for triangles of varying orders. These are not tensor product of 
    /// simple <see cref="GaussLegendre1D_old"/> rules. Instead the points are derived by symmetric Gaussian quadratures for
    /// triangles which are more efficient. 
    /// The theory is published in: "High degree efficient symmetrical Gaussian quadrature rules for the triangle", 
    /// D. A. Dunavant, 1985. It is also explained in http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF. The
    /// integration rules listed in the book: "Finite Element Procedures", K. J. Bathe, 2014 are a subset of the Dunavant rules
    /// (Gauss points are ordered differently), even though Bathe attributes them to G. R. Cowper.
    /// Note that in the integral: \iint_T f(xi, eta)*dxi*deta = \sum_{i=1}^N w_i * f(xi_i, eta_i), the factor 1/2 that is
    /// commonly seen in front of the sum in some sources, is incorporated into the weight factor here.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public sealed class TriangleQuadratureSymmetricGaussian : IQuadrature2D
    {
        public static readonly TriangleQuadratureSymmetricGaussian Order1Point1 = new TriangleQuadratureSymmetricGaussian(
            new GaussPoint2D(0.33333333333333, 0.33333333333333, 0.5000000000000000));

        // WARNING: Do not change their order. ExtrapolationGaussTriangular3Points depends on it.
        public static readonly TriangleQuadratureSymmetricGaussian Order2Points3 = new TriangleQuadratureSymmetricGaussian(
            new GaussPoint2D(0.66666666666667, 0.16666666666667, 0.5 * 0.33333333333333),
            new GaussPoint2D(0.16666666666667, 0.66666666666667, 0.5 * 0.33333333333333),
            new GaussPoint2D(0.16666666666667, 0.16666666666667, 0.5 * 0.33333333333333));

        public static readonly TriangleQuadratureSymmetricGaussian Order3Points4 = new TriangleQuadratureSymmetricGaussian(
            new GaussPoint2D(0.33333333333333, 0.33333333333333, 0.5 * -0.5625000000000),
            new GaussPoint2D(0.20000000000000, 0.20000000000000, 0.5 * 0.52083333333333),
            new GaussPoint2D(0.20000000000000, 0.60000000000000, 0.5 * 0.52083333333333),
            new GaussPoint2D(0.60000000000000, 0.20000000000000, 0.5 * 0.52083333333333));

        public static readonly TriangleQuadratureSymmetricGaussian Order4Points6 = new TriangleQuadratureSymmetricGaussian(
            new GaussPoint2D(0.44594849091597, 0.44594849091597, 0.5 * 0.22338158967801),
            new GaussPoint2D(0.44594849091597, 0.10810301816807, 0.5 * 0.22338158967801),
            new GaussPoint2D(0.10810301816807, 0.44594849091597, 0.5 * 0.22338158967801),
            new GaussPoint2D(0.09157621350977, 0.09157621350977, 0.5 * 0.10995174365532),
            new GaussPoint2D(0.09157621350977, 0.81684757298046, 0.5 * 0.10995174365532),
            new GaussPoint2D(0.81684757298046, 0.09157621350977, 0.5 * 0.10995174365532));

        public static readonly TriangleQuadratureSymmetricGaussian Order5Points7 = new TriangleQuadratureSymmetricGaussian(
            new GaussPoint2D(0.33333333333333, 0.33333333333333, 0.5 * 0.22500000000000),
            new GaussPoint2D(0.47014206410511, 0.47014206410511, 0.5 * 0.13239415278851),
            new GaussPoint2D(0.47014206410511, 0.05971587178977, 0.5 * 0.13239415278851),
            new GaussPoint2D(0.05971587178977, 0.47014206410511, 0.5 * 0.13239415278851),
            new GaussPoint2D(0.10128650732346, 0.10128650732346, 0.5 * 0.12593918054483),
            new GaussPoint2D(0.10128650732346, 0.79742698535309, 0.5 * 0.12593918054483),
            new GaussPoint2D(0.79742698535309, 0.10128650732346, 0.5 * 0.12593918054483));

        public static readonly TriangleQuadratureSymmetricGaussian Order6Points12 = new TriangleQuadratureSymmetricGaussian(
            new GaussPoint2D(0.24928674517091, 0.24928674517091, 0.5 * 0.11678627572638),
            new GaussPoint2D(0.24928674517091, 0.50142650965818, 0.5 * 0.11678627572638),
            new GaussPoint2D(0.50142650965818, 0.24928674517091, 0.5 * 0.11678627572638),
            new GaussPoint2D(0.06308901449150, 0.06308901449150, 0.5 * 0.05084490637021),
            new GaussPoint2D(0.06308901449150, 0.87382197101700, 0.5 * 0.05084490637021),
            new GaussPoint2D(0.87382197101700, 0.06308901449150, 0.5 * 0.05084490637021),
            new GaussPoint2D(0.31035245103378, 0.63650249912140, 0.5 * 0.08285107561837),
            new GaussPoint2D(0.63650249912140, 0.05314504984482, 0.5 * 0.08285107561837),
            new GaussPoint2D(0.05314504984482, 0.31035245103378, 0.5 * 0.08285107561837),
            new GaussPoint2D(0.63650249912140, 0.31035245103378, 0.5 * 0.08285107561837),
            new GaussPoint2D(0.31035245103378, 0.05314504984482, 0.5 * 0.08285107561837),
            new GaussPoint2D(0.05314504984482, 0.63650249912140, 0.5 * 0.08285107561837));

        public static readonly TriangleQuadratureSymmetricGaussian Order7Points13 = new TriangleQuadratureSymmetricGaussian(
            new GaussPoint2D(0.33333333333333, 0.33333333333333, 0.5 * -0.14957004446768),
            new GaussPoint2D(0.26034596607904, 0.26034596607904, 0.5 * 0.17561525743321),
            new GaussPoint2D(0.26034596607904, 0.47930806784192, 0.5 * 0.17561525743321),
            new GaussPoint2D(0.47930806784192, 0.26034596607904, 0.5 * 0.17561525743321),
            new GaussPoint2D(0.06513010290222, 0.06513010290222, 0.5 * 0.05334723560884),
            new GaussPoint2D(0.06513010290222, 0.86973979419557, 0.5 * 0.05334723560884),
            new GaussPoint2D(0.86973979419557, 0.06513010290222, 0.5 * 0.05334723560884),
            new GaussPoint2D(0.31286549600487, 0.63844418856981, 0.5 * 0.07711376089026),
            new GaussPoint2D(0.63844418856981, 0.04869031542532, 0.5 * 0.07711376089026),
            new GaussPoint2D(0.04869031542532, 0.31286549600487, 0.5 * 0.07711376089026),
            new GaussPoint2D(0.63844418856981, 0.31286549600487, 0.5 * 0.07711376089026),
            new GaussPoint2D(0.31286549600487, 0.04869031542532, 0.5 * 0.07711376089026),
            new GaussPoint2D(0.04869031542532, 0.63844418856981, 0.5 * 0.07711376089026));

        public static readonly TriangleQuadratureSymmetricGaussian Order8Points16 = new TriangleQuadratureSymmetricGaussian(
            new GaussPoint2D(0.33333333333333, 0.33333333333333, 0.5 * 0.14431560767779),
            new GaussPoint2D(0.45929258829272, 0.45929258829272, 0.5 * 0.09509163426728),
            new GaussPoint2D(0.45929258829272, 0.08141482341455, 0.5 * 0.09509163426728),
            new GaussPoint2D(0.08141482341455, 0.45929258829272, 0.5 * 0.09509163426728),
            new GaussPoint2D(0.17056930775176, 0.17056930775176, 0.5 * 0.10321737053472),
            new GaussPoint2D(0.17056930775176, 0.65886138449648, 0.5 * 0.10321737053472),
            new GaussPoint2D(0.65886138449648, 0.17056930775176, 0.5 * 0.10321737053472),
            new GaussPoint2D(0.05054722831703, 0.05054722831703, 0.5 * 0.03245849762320),
            new GaussPoint2D(0.05054722831703, 0.89890554336594, 0.5 * 0.03245849762320),
            new GaussPoint2D(0.89890554336594, 0.05054722831703, 0.5 * 0.03245849762320),
            new GaussPoint2D(0.26311282963464, 0.72849239295540, 0.5 * 0.02723031417443),
            new GaussPoint2D(0.72849239295540, 0.00839477740996, 0.5 * 0.02723031417443),
            new GaussPoint2D(0.00839477740996, 0.26311282963464, 0.5 * 0.02723031417443),
            new GaussPoint2D(0.72849239295540, 0.26311282963464, 0.5 * 0.02723031417443),
            new GaussPoint2D(0.26311282963464, 0.00839477740996, 0.5 * 0.02723031417443),
            new GaussPoint2D(0.00839477740996, 0.72849239295540, 0.5 * 0.02723031417443));

        private TriangleQuadratureSymmetricGaussian(params GaussPoint2D[] points)
        {
            this.IntegrationPoints = new List<GaussPoint2D>(points);
        }

        public IReadOnlyList<GaussPoint2D> IntegrationPoints { get; }
    }
}
