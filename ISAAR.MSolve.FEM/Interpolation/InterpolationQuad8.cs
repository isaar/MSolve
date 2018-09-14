using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;

// Quad8 nodes:
// 3 -- 6 -- 2
// |         |
// 7         5
// |         |
// 0 -- 4 -- 1

namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Isoparametric interpolation of a quadrilateral finite element with 8 nodes. Quadratic shape functions. 
    /// Implements Singleton pattern.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class InterpolationQuad8 : IsoparametricInterpolation2DBase
    {
        private static readonly InterpolationQuad8 uniqueInstance = new InterpolationQuad8();

        private InterpolationQuad8() : base(CellType2D.Quad8, 8)
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
            };
        }

        /// <summary>
        /// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
        /// nodes matches the order of the shape functions and is always the same for each element.
        /// </summary>
        public override IReadOnlyList<NaturalPoint2D> NodalNaturalCoordinates { get; }

        /// <summary>
        /// Get the unique <see cref="InterpolationQuad8"/> object for the whole program. Thread safe.
        /// </summary>
        public static InterpolationQuad8 UniqueInstance => uniqueInstance;

        /// <summary>
        /// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <returns></returns>
        public override IInverseInterpolation2D CreateInverseMappingFor(IReadOnlyList<Node2D> nodes)
            => throw new NotImplementedException("Requires an iterative procedure.");

        protected override sealed double[] EvaluateAt(double xi, double eta)
        {
            //TODO: cache some quantities, e.g. 1-xi, 1+xi, etc.
            double xi_2 = xi * xi;
            double eta_2 = eta * eta;

            var values = new double[8];
            values[0] = 0.25 * (1 - xi) * (1 - eta) * (-xi - eta - 1);
            values[1] = 0.25 * (1 + xi) * (1 - eta) * (xi - eta - 1);
            values[2] = 0.25 * (1 + xi) * (1 + eta) * (xi + eta - 1);
            values[3] = 0.25 * (1 - xi) * (1 + eta) * (-xi + eta - 1);

            values[4] = 0.5 * (1 - xi_2) * (1 - eta);
            values[5] = 0.5 * (1 + xi) * (1 - eta_2);
            values[6] = 0.5 * (1 - xi_2) * (1 + eta);
            values[7] = 0.5 * (1 - xi) * (1 - eta_2);

            return values;
        }

        protected override sealed double[,] EvaluateGradientsAt(double xi, double eta)
        {
            //TODO: cache some quantities, e.g. 1-xi, 1+xi, etc.
            double xi2 = xi * 2.0;
            double eta2 = eta * 2.0;
            double xiSq = xi * xi;
            double etaSq = eta * eta;
            double xiEta = xi * eta;

            var derivatives = new double[8, 2];
            derivatives[0, 0] = -0.25 * (xi2 + eta) * (eta - 1);
            derivatives[0, 1] = -0.25 * (xi - 1) * (xi + eta2);
            derivatives[1, 0] = -0.25 * (xi2 - eta) * (eta - 1);
            derivatives[1, 1] = -0.25 * (xi + 1) * (xi - eta2);
            derivatives[2, 0] = -0.25 * (xi2 + eta) * (-eta - 1);
            derivatives[2, 1] = -0.25 * (xi + 1) * (-xi - eta2);
            derivatives[3, 0] = -0.25 * (xi2 - eta) * (-eta - 1);
            derivatives[3, 1] = -0.25 * (xi - 1) * (-xi + eta2);

            derivatives[4, 0] = xi * (eta - 1);
            derivatives[4, 1] = 0.5 * (xiSq - 1);
            derivatives[5, 0] = 0.5 * (1 - etaSq);
            derivatives[5, 1] = -eta * (xi + 1);
            derivatives[6, 0] = -xi * (eta + 1);
            derivatives[6, 1] = -0.5 * (xiSq - 1);
            derivatives[7, 0] = -0.5 * (1 - etaSq);
            derivatives[7, 1] = eta * (xi - 1);
            return derivatives;
        }
    }
}
