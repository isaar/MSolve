using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Geometry.Coordinates
{
    /// <summary>
    /// Point in a 2-dimensional cartesian coordinate system, which is local to a finite element. Immutable.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class NaturalPoint2D : IPoint2D
    {
        protected readonly double xi;
        protected readonly double eta;

        /// <summary>
        /// Instantiates a <see cref="NaturalPoint2D"/>.
        /// </summary>
        /// <param name="xi">The coordinate of the point along local axis Xi.</param>
        /// <param name="eta">The coordinate of the point along local axis Eta.</param>
        public NaturalPoint2D(double xi, double eta)
        {
            this.xi = xi;
            this.eta = eta;
        }

        /// <summary>
        /// Instantiates a <see cref="NaturalPoint2D"/>.
        /// </summary>
        /// <param name="coordinates">Vector with the coordinates of the point. Length = 2.</param>
        public NaturalPoint2D(double[] coordinates)
        {
            this.xi = coordinates[0];
            this.eta = coordinates[1];
        }

        /// <summary>
        /// Vector with the coordinates of the point. Length = 2.
        /// </summary>
        public double[] Coordinates => new double[] { xi, eta };

        /// <summary>
        /// The coordinate of the point along local axis Xi.
        /// </summary>
        public double Xi => xi;

        /// <summary>
        /// The coordinate of the point along local axis Eta.
        /// </summary>
        public double Eta => eta;

        public double X1 => xi;

        public double X2 => eta;

        public override string ToString() => $"(xi, eta) = ({xi}, {eta})";
    }
}
