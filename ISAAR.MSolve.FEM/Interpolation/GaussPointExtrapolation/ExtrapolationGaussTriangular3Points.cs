using ISAAR.MSolve.Geometry.Coordinates;
using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Quadratures;

namespace ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation
{
    /// <summary>
    /// Calculates extrapolations of scalar, vector and tensor fields from the integration points of 
    /// <see cref="GaussQuadratureForTriangles.Order2Points3"/>. This can be done for any point, but utility methods for directly 
    /// outputting the extrapolated fields at the nodes of finite elements are also provided.
    /// Implements Singleton pattern.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ExtrapolationGaussTriangular3Points: GaussPointExtrapolation2DBase
    {
        private const double oneOverThree = 1.0 / 3.0;
        private static readonly ExtrapolationGaussTriangular3Points uniqueInstance = new ExtrapolationGaussTriangular3Points();

        private ExtrapolationGaussTriangular3Points() : base(GaussQuadratureForTriangles.Order2Points3)
        { }

        /// <summary>
        /// Get the unique <see cref="ExtrapolationGaussTriangular3Points"/> object for the whole program. Thread safe.
        /// </summary>
        public static ExtrapolationGaussTriangular3Points UniqueInstance => uniqueInstance;

        protected override double[] EvaluateExtrapolationFunctionsAt(NaturalPoint2D point)
        {
            // Coordinates of the point in the auxiliary coordinate system of an imaginary "Gauss element" that has the Gauss 
            // points as its nodes.
            double r = 2.0 * point.Xi - oneOverThree;
            double s = 2.0 * point.Eta - oneOverThree;

            // Shape functions of the imaginary "Gauss element".
            // Each shape function corresponds to an integration point of GaussQuadratureForTriangles.Order2Points3. Therefore 
            // their order must be the same: point at right angle, point on Xi, point on Eta.
            var shapeFunctions = new double[3];
            shapeFunctions[0] = 1 - r - s;
            shapeFunctions[1] = r;
            shapeFunctions[2] = s;
            return shapeFunctions;
        }
    }
}
