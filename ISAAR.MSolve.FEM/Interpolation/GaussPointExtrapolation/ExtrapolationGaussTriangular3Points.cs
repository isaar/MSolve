using ISAAR.MSolve.FEM.Integration.Quadratures;
using ISAAR.MSolve.Geometry.Coordinates;
using System;
using System.Collections.Generic;
using System.Text;

// TODO: This order of shape functions produces smooth plots and agrees with Abaqus. But why?
// This extrapolation is used for Tri6 interpolation with 3 Gauss points, where the numbering is:
//
// eta
// ^
// | s
// | ^
//   |
// 1
// | \
// | 1 \
// |     \
// 4       3
// |         \
// | 2      0  \    --> r
// |             \
// 2 ---- 5 ----- 0   --> xi
//
// Shouldn't the shape functions be N0 = r, N1 = s, N2 = 1 - r - s? When I tried it, the tensors at nodes where 
// permuted with respect to Abaqus and the plots (after averaging) were not smooth (e.g. for a 2D cantilever under 
// bending).
// Note that the Gauss point order in Abaqus is identical to GaussQuadratureForTrianglesSymmetric.Order2Points3:
// A0 = (2/3, 1/6), A1 = (1/6, 2/3), A2 = (1/6, 1/6).
// Perhaps I should evaluate the shape functions at each Gauss point to see what is happening.
namespace ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation
{
    /// <summary>
    /// Calculates extrapolations of scalar, vector and tensor fields from the integration points of 
    /// <see cref="GaussQuadratureForTrianglesSymmetric.Order2Points3"/>. This can be done for any point, but utility methods for directly 
    /// outputting the extrapolated fields at the nodes of finite elements are also provided.
    /// Implements Singleton pattern.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ExtrapolationGaussTriangular3Points: GaussPointExtrapolation2DBase
    {
        private const double oneOverThree = 1.0 / 3.0;
        private static readonly ExtrapolationGaussTriangular3Points uniqueInstance = new ExtrapolationGaussTriangular3Points();

        private ExtrapolationGaussTriangular3Points() : base(GaussQuadratureForTrianglesSymmetric.Order2Points3)
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
            // Each shape function corresponds to an integration point of GaussQuadratureForTrianglesSymmetric.Order2Points3.  
            // Therefore their order must be the same: point on Xi, point on Eta, point at right angle. This might differ from 
            // the node order in InterpolationTri3.
            var shapeFunctions = new double[3];
            shapeFunctions[0] = 1 - r - s;
            shapeFunctions[1] = r;
            shapeFunctions[2] = s;
            return shapeFunctions;
        }
    }
}
