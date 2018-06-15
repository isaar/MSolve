using ISAAR.MSolve.FEM.Integration.Quadratures;
using ISAAR.MSolve.Geometry.Coordinates;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation
{
    /// <summary>
    /// Calculates extrapolations of scalar, vector and tensor fields from the integration points of 2-by-2 Gauss-Legendre 
    /// quadrature. This can be done for any point, but utility methods for directly outputting the extrapolated fields at the
    /// nodes of finite elements are also provided.
    /// Implements Singleton pattern.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ExtrapolationGaussLegendre2x2: GaussPointExtrapolation2DBase
    {
        private static readonly double sqrt3 = Math.Sqrt(3.0);
        private static readonly ExtrapolationGaussLegendre2x2 uniqueInstance = new ExtrapolationGaussLegendre2x2();

        private ExtrapolationGaussLegendre2x2(): base(GaussLegendre2D.Order2x2)
        { }

        /// <summary>
        /// Get the unique <see cref="ExtrapolationGaussLegendre2x2"/> object for the whole program. Thread safe.
        /// </summary>
        public static ExtrapolationGaussLegendre2x2 UniqueInstance { get { return uniqueInstance; } }

        protected override double[] EvaluateExtrapolationFunctionsAt(NaturalPoint2D point)
        {
            // Coordinates of the point in the auxiliary coordinate system of an imaginary "Gauss element" that has the Gauss 
            // points as its nodes.
            double r = sqrt3 * point.Xi;
            double s = sqrt3 * point.Eta;

            // Shape functions of the imaginary "Gauss element".
            var shapeFunctions = new double[4];
            shapeFunctions[0] = 0.25 * (1 - r) * (1 - s);
            shapeFunctions[1] = 0.25 * (1 + r) * (1 - s);
            shapeFunctions[2] = 0.25 * (1 + r) * (1 + s);
            shapeFunctions[3] = 0.25 * (1 - r) * (1 + s);
            return shapeFunctions;
        }
    }
}
