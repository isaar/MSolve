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
    public class InterpolationQuad4: IsoparametricInterpolation2DBase
    {
        private static readonly InterpolationQuad4 uniqueInstance = new InterpolationQuad4();

        private InterpolationQuad4(): base(CellType2D.Quad4, 4)
        {
            NodalNaturalCoordinates = new NaturalPoint2D[]
            {
                new NaturalPoint2D(-1.0, -1.0),
                new NaturalPoint2D(+1.0, -1.0),
                new NaturalPoint2D(+1.0, +1.0),
                new NaturalPoint2D(-1.0, +1.0)
            };
        }

        /// <summary>
        /// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
        /// nodes matches the order of the shape functions and is always the same for each element.
        /// </summary>
        public override IReadOnlyList<NaturalPoint2D> NodalNaturalCoordinates { get; }

        /// <summary>
        /// Get the unique <see cref="InterpolationQuad4"/> object for the whole program. Thread safe.
        /// </summary>
        public static InterpolationQuad4 UniqueInstance => uniqueInstance;

        /// <summary>
        /// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <returns></returns>
        public override IInverseInterpolation2D CreateInverseMappingFor(IReadOnlyList<Node2D> nodes)
            => new InverseInterpolationQuad4(nodes);

        protected override sealed double[] EvaluateAt(double xi, double eta)
        {
            var values = new double[4];
            values[0] = 0.25 * (1 - xi) * (1 - eta);
            values[1] = 0.25 * (1 + xi) * (1 - eta);
            values[2] = 0.25 * (1 + xi) * (1 + eta);
            values[3] = 0.25 * (1 - xi) * (1 + eta);
            return values;
        }

        protected override sealed double[,] EvaluateGradientsAt(double xi, double eta)
        {
            var derivatives = new double[4, 2];
            derivatives[0, 0] = -0.25 * (1 - eta);
            derivatives[0, 1] = -0.25 * (1 - xi);
            derivatives[1, 0] = 0.25 * (1 - eta);
            derivatives[1, 1] = -0.25 * (1 + xi);
            derivatives[2, 0] = 0.25 * (1 + eta);
            derivatives[2, 1] = 0.25 * (1 + xi);
            derivatives[3, 0] = -0.25 * (1 + eta);
            derivatives[3, 1] = 0.25 * (1 - xi);
            return derivatives;
        }
    }
}
