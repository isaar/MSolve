using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Points;

namespace ISAAR.MSolve.Discretization.Integration.Quadratures
{
    /// <summary>
    /// Enum class with the 1D Gauss-Legendre integration rules of varying orders. 
    /// For more see https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Legendre_quadrature.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public sealed class GaussLegendre1D : IQuadrature1D
    {
        public static readonly GaussLegendre1D Order1 = new GaussLegendre1D(
            new GaussPoint1D(0.0, 2.0));

        public static readonly GaussLegendre1D Order2 = new GaussLegendre1D(
            new GaussPoint1D(-0.5773502691896257645091488, 1.0),
            new GaussPoint1D(+0.5773502691896257645091488, 1.0));

        public static readonly GaussLegendre1D Order3 = new GaussLegendre1D(
            new GaussPoint1D(-0.7745966692414833770358531, 0.5555555555555555555555556),
            new GaussPoint1D(0.0, 0.8888888888888888888888889),
            new GaussPoint1D(+0.7745966692414833770358531, 0.5555555555555555555555556));

        public static readonly GaussLegendre1D Order4 = new GaussLegendre1D(
            new GaussPoint1D(-0.8611363115940525752239465, 0.3478548451374538573730639),
            new GaussPoint1D(-0.3399810435848562648026658, 0.6521451548625461426269361),
            new GaussPoint1D(+0.3399810435848562648026658, 0.6521451548625461426269361),
            new GaussPoint1D(+0.8611363115940525752239465, 0.3478548451374538573730639));

        public static readonly GaussLegendre1D Order5 = new GaussLegendre1D(
            new GaussPoint1D(-0.9061798459386639927976269, 0.2369268850561890875142640),
            new GaussPoint1D(-0.5384693101056830910363144, 0.4786286704993664680412915),
            new GaussPoint1D(0, 0.5688888888888888888888889),
            new GaussPoint1D(+0.5384693101056830910363144, 0.4786286704993664680412915),
            new GaussPoint1D(+0.9061798459386639927976269, 0.2369268850561890875142640));

        public static readonly GaussLegendre1D Order6 = new GaussLegendre1D(
            new GaussPoint1D(-0.9324695142031520278123016, 0.1713244923791703450402961),
            new GaussPoint1D(-0.6612093864662645136613996, 0.3607615730481386075698335),
            new GaussPoint1D(-0.2386191860831969086305017, 0.4679139345726910473898703),
            new GaussPoint1D(+0.2386191860831969086305017, 0.4679139345726910473898703),
            new GaussPoint1D(+0.6612093864662645136613996, 0.3607615730481386075698335),
            new GaussPoint1D(+0.9324695142031520278123016, 0.1713244923791703450402961));

        public static readonly GaussLegendre1D Order7 = new GaussLegendre1D(
            new GaussPoint1D(-0.9491079123427585245261897, 0.1294849661688696932706114),
            new GaussPoint1D(-0.7415311855993944398638648, 0.2797053914892766679014678),
            new GaussPoint1D(-0.4058451513773971669066064, 0.3818300505051189449503698),
            new GaussPoint1D(0.0, 0.4179591836734693877551020),
            new GaussPoint1D(+0.4058451513773971669066064, 0.3818300505051189449503698),
            new GaussPoint1D(+0.7415311855993944398638648, 0.2797053914892766679014678),
            new GaussPoint1D(+0.9491079123427585245261897, 0.1294849661688696932706114));

        public static readonly GaussLegendre1D Order8 = new GaussLegendre1D(
            new GaussPoint1D(-0.9602898564975362316835609, 0.1012285362903762591525314),
            new GaussPoint1D(-0.7966664774136267395915539, 0.2223810344533744705443560),
            new GaussPoint1D(-0.5255324099163289858177390, 0.3137066458778872873379622),
            new GaussPoint1D(-0.1834346424956498049394761, 0.3626837833783619829651504),
            new GaussPoint1D(+0.1834346424956498049394761, 0.3626837833783619829651504),
            new GaussPoint1D(+0.5255324099163289858177390, 0.3137066458778872873379622),
            new GaussPoint1D(+0.7966664774136267395915539, 0.2223810344533744705443560),
            new GaussPoint1D(+0.9602898564975362316835609, 0.1012285362903762591525314));

        public static readonly GaussLegendre1D Order9 = new GaussLegendre1D(
            new GaussPoint1D(-0.9681602395076260898355762, 0.0812743883615744119718922),
            new GaussPoint1D(-0.8360311073266357942994298, 0.1806481606948574040584720),
            new GaussPoint1D(-0.6133714327005903973087020, 0.2606106964029354623187429),
            new GaussPoint1D(-0.3242534234038089290385380, 0.3123470770400028400686304),
            new GaussPoint1D(0, 0.3302393550012597631645251),
            new GaussPoint1D(+0.3242534234038089290385380, 0.3123470770400028400686304),
            new GaussPoint1D(+0.6133714327005903973087020, 0.2606106964029354623187429),
            new GaussPoint1D(+0.8360311073266357942994298, 0.1806481606948574040584720),
            new GaussPoint1D(+0.9681602395076260898355762, 0.0812743883615744119718922));

        public static readonly GaussLegendre1D Order10 = new GaussLegendre1D(
            new GaussPoint1D(-0.9739065285171717200779640, 0.0666713443086881375935688),
            new GaussPoint1D(-0.8650633666889845107320967, 0.1494513491505805931457763),
            new GaussPoint1D(-0.6794095682990244062343274, 0.2190863625159820439955349),
            new GaussPoint1D(-0.1488743389816312108848260, 0.2955242247147528701738930),
            new GaussPoint1D(-0.4333953941292471907992659, 0.2692667193099963550912269),
            new GaussPoint1D(+0.1488743389816312108848260, 0.2955242247147528701738930),
            new GaussPoint1D(+0.4333953941292471907992659, 0.2692667193099963550912269),
            new GaussPoint1D(+0.6794095682990244062343274, 0.2190863625159820439955349),
            new GaussPoint1D(+0.8650633666889845107320967, 0.1494513491505805931457763),
            new GaussPoint1D(+0.9739065285171717200779640, 0.0666713443086881375935688));

        private GaussLegendre1D(params GaussPoint1D[] points)
        {
            IntegrationPoints = new List<GaussPoint1D>(points);
        }

        /// <summary>
        /// The integrations points are sorted in increasing xi order. This order is strictly defined for each quadrature and 
        /// cannot change.
        /// </summary>
        public IReadOnlyList<GaussPoint1D> IntegrationPoints { get; }
    }
}
