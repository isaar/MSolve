using System;
using System.Collections.Generic;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.LinearAlgebra.Matrices;

// Quad9 nodes:
// 3 -- 6 -- 2
// |    |    |
// 7 -- 8 -- 5
// |    |    |
// 0 -- 4 -- 1

namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Isoparametric interpolation of a quadrilateral finite element with 9 nodes. Quadratic shape functions. 
    /// Implements Singleton pattern.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class InterpolationQuad9 : IsoparametricInterpolation2DBase
    {
        private static readonly InterpolationQuad9 uniqueInstance = new InterpolationQuad9();

        private InterpolationQuad9() : base(CellType.Quad9, 9)
        {
            NodalNaturalCoordinates = new NaturalPoint2D[]
            {
                new NaturalPoint2D(-1.0, -1.0),
                new NaturalPoint2D(+1.0, -1.0),
                new NaturalPoint2D(+1.0, +1.0),
                new NaturalPoint2D(-1.0, +1.0),

                new NaturalPoint2D(+0.0, -1.0),
                new NaturalPoint2D(+1.0, +0.0),
                new NaturalPoint2D(+0.0, +1.0),
                new NaturalPoint2D(-1.0, +0.0),

                new NaturalPoint2D(+0.0, +0.0)
            };
        }

        /// <summary>
        /// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
        /// nodes matches the order of the shape functions and is always the same for each element.
        /// </summary>
        public override IReadOnlyList<NaturalPoint2D> NodalNaturalCoordinates { get; }

        /// <summary>
        /// Get the unique <see cref="InterpolationQuad9"/> object for the whole program. Thread safe.
        /// </summary>
        public static InterpolationQuad9 UniqueInstance => uniqueInstance;

        /// <summary>
        /// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <returns></returns>
        public override IInverseInterpolation2D CreateInverseMappingFor(IReadOnlyList<Node_v2> nodes)
            => throw new NotImplementedException("Requires an iterative procedure.");

        protected override sealed double[] EvaluateAt(double xi, double eta)
        {
            double xiEtaOver4 = 0.25 * xi * eta;
            double xi_2 = xi * xi;
            double eta_2 = eta * eta;

            var values = new double[9];
            values[0] = xiEtaOver4 * (1 - xi) * (1 - eta);
            values[1] = -xiEtaOver4 * (1 + xi) * (1 - eta);
            values[2] = xiEtaOver4 * (1 + xi) * (1 + eta);
            values[3] = -xiEtaOver4 * (1 - xi) * (1 + eta);

            values[4] = -0.5 * eta * (1 - xi_2) * (1 - eta);
            values[5] = 0.5 * xi * (1 + xi) * (1 - eta_2);
            values[6] = 0.5 * eta * (1 - xi_2) * (1 + eta);
            values[7] = -0.5 * xi * (1 - xi) * (1 - eta_2);

            values[8] = (1 - xi_2) * (1 - eta_2);
            return values;
        }

        protected override sealed Matrix EvaluateGradientsAt(double xi, double eta)
        {
            double xi2 = xi * 2.0;
            double eta2 = eta * 2.0;
            double xiSq = xi * xi;
            double etaSq = eta * eta;
            double xiEta = xi * eta;

            var derivatives = Matrix.CreateZero(9, 2);
            derivatives[0, 0] = 0.25 * eta * (1 - xi2) * (1 - eta);
            derivatives[0, 1] = 0.25 * xi * (1 - xi) * (1 - eta2);
            derivatives[1, 0] = -0.25 * eta * (1 + xi2) * (1 - eta);
            derivatives[1, 1] = -0.25 * xi * (1 + xi) * (1 - eta2);
            derivatives[2, 0] = 0.25 * eta * (1 + xi2) * (1 + eta);
            derivatives[2, 1] = 0.25 * xi * (1 + xi) * (1 + eta2);
            derivatives[3, 0] = -0.25 * eta * (1 - xi2) * (1 + eta);
            derivatives[3, 1] = -0.25 * xi * (1 - xi) * (1 + eta2);

            derivatives[4, 0] = xiEta * (1 - eta);
            derivatives[4, 1] = -0.5 * (1 - xiSq) * (1 - eta2);
            derivatives[5, 0] = 0.5 * (1 + xi2) * (1 - etaSq);
            derivatives[5, 1] = -xiEta * (1 + xi);
            derivatives[6, 0] = -xiEta * (1 + eta);
            derivatives[6, 1] = 0.5 * (1 - xiSq) * (1 + eta2);
            derivatives[7, 0] = -0.5 * (1 - xi2) * (1 - etaSq);
            derivatives[7, 1] = xiEta * (1 - xi);

            derivatives[8, 0] = -2 * xi * (1 - etaSq);
            derivatives[8, 1] = -2 * eta * (1 - xiSq);
            return derivatives;
        }
    }
}
