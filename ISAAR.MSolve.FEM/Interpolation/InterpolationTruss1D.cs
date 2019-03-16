using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;

// Quad4 nodes:
// 3 -- 2
// |    |
// 0 -- 1

namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Isoparametric interpolation of a quadrilateral finite element with 4 nodes. Linear shape functions.
    /// Implements Singleton pattern.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class InterpolationTruss1D : IsoparametricInterpolation2DBase
    {
        private static readonly InterpolationTruss1D uniqueInstance = new InterpolationTruss1D();

        private InterpolationTruss1D(): base(CellType2D.Quad4, 4)
        {
            NodalNaturalCoordinates = new NaturalPoint2D[]
            {
                new NaturalPoint2D(-1.0, 0.0),
                new NaturalPoint2D(+1.0, 0.0)
            };
        }

        /// <summary>
        /// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
        /// nodes matches the order of the shape functions and is always the same for each element.
        /// </summary>
        public override IReadOnlyList<NaturalPoint2D> NodalNaturalCoordinates { get; }

        /// <summary>
        /// Get the unique <see cref="InterpolationTruss1D"/> object for the whole program. Thread safe.
        /// </summary>
        public static InterpolationTruss1D UniqueInstance => uniqueInstance;

        /// <summary>
        /// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <returns></returns>
        public override IInverseInterpolation2D CreateInverseMappingFor(IReadOnlyList<Node2D> nodes)
            => new InverseInterpolationTruss1D(nodes);

        protected override sealed double[] EvaluateAt(double xi, double eta)
        {
            var values = new double[4];
            values[0] = 0.50 * (1 - xi);
            values[1] = 0.50 * (1 + xi);
            return values;
        }

        protected override sealed double[,] EvaluateGradientsAt(double xi, double eta)
        {
            var derivatives = new double[1, 2];
            derivatives[0, 0] = -0.50; // N1,ksi
            derivatives[0, 1] = +0.50; // N2,ksi
            return derivatives;
        }
    }
}
