using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;

// Truss nodes:
// 0 -- 1

namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Isoparametric interpolation of a quadrilateral finite element with 4 nodes. Linear shape functions.
    /// Implements Singleton pattern.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class InterpolationTruss1D
    {
        private static readonly InterpolationTruss1D uniqueInstance = new InterpolationTruss1D();

        private InterpolationTruss1D()
        {
            NodalNaturalCoordinates = new NaturalPoint1D[]
            {
                new NaturalPoint1D(-1.0),
                new NaturalPoint1D(+1.0)
            };
        }

        /// <summary>
        /// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
        /// nodes matches the order of the shape functions and is always the same for each element.
        /// </summary>
        public IReadOnlyList<NaturalPoint1D> NodalNaturalCoordinates { get; }

        /// <summary>
        /// Get the unique <see cref="InterpolationTruss1D"/> object for the whole program. Thread safe.
        /// </summary>
        public static InterpolationTruss1D UniqueInstance => uniqueInstance;

        /// <summary>
        /// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <returns></returns>
        //public override IInverseInterpolation1D CreateInverseMappingFor(IReadOnlyList<Node> nodes)
        //    => new InverseInterpolationTruss1D(nodes);

        protected double[] EvaluateAt(double xi)
        {
            var values = new double[2];
            values[0] = 0.50 * (1 - xi);
            values[1] = 0.50 * (1 + xi);
            return values;
        }

        protected double[,] EvaluateGradientsAt(double xi)
        {
            var derivatives = new double[1, 2];
            derivatives[0, 0] = -0.50; // N1,ksi
            derivatives[0, 1] = +0.50; // N2,ksi
            return derivatives;
        }
    }
}
