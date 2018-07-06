using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;

// Tri3 nodes:
// 1
// |  \
// 2 -- 0

//TODO: See https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch24.d/IFEM.Ch24.pdf for optimizations
namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Isoparametric interpolation of a triangular finite element with 3 nodes. Linear shape functions. 
    /// Implements Singleton pattern.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class InterpolationTri3 : IsoparametricInterpolation2DBase
    {
        private static readonly InterpolationTri3 uniqueInstance = new InterpolationTri3();

        private InterpolationTri3() : base(CellType2D.Tri3, 3)
        {
            NodalNaturalCoordinates = new NaturalPoint2D[]
            {
                new NaturalPoint2D(1.0, 0.0),
                new NaturalPoint2D(0.0, 1.0),
                new NaturalPoint2D(0.0, 0.0)
            };
        }

        /// <summary>
        /// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
        /// nodes matches the order of the shape functions and is always the same for each element.
        /// </summary>
        public override IReadOnlyList<NaturalPoint2D> NodalNaturalCoordinates { get; }

        /// <summary>
        /// Get the unique <see cref="InterpolationTri3"/> object for the whole program. Thread safe.
        /// </summary>
        public static InterpolationTri3 UniqueInstance => uniqueInstance;

        /// <summary>
        /// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <returns></returns>
        public override IInverseInterpolation2D CreateInverseMappingFor(IReadOnlyList<Node2D> nodes)
            => new InverseInterpolationTri3(nodes);
        
        protected override sealed double[] EvaluateAt(double xi, double eta)
        {
            var values = new double[3];
            values[0] = xi;
            values[1] = eta;
            values[2] = 1 - xi - eta;
            return values;
        }

        protected override sealed double[,] EvaluateGradientsAt(double xi, double eta)
        {
            var derivatives = new double[3, 2];
            derivatives[0, 0] = +1.0;
            derivatives[0, 1] = +0.0;
            derivatives[1, 0] = +0.0;
            derivatives[1, 1] = +1.0;
            derivatives[2, 0] = -1.0;
            derivatives[2, 1] = -1.0;
            return derivatives;
        }
    }
}
