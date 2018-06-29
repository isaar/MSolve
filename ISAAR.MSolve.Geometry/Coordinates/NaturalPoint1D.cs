using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Geometry.Coordinates
{
    /// <summary>
    /// Point in a 1-dimensional cartesian coordinate system, which is local to a finite element. Immutable.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class NaturalPoint1D : IPoint1D
    {
        protected readonly double xi;

        /// <summary>
        /// Instantiates a <see cref="NaturalPoint1D"/>.
        /// </summary>
        /// <param name="xi">The coordinate of the point along the single axis Xi.</param>
        public NaturalPoint1D(double xi)
        {
            this.xi = xi;
        }

        /// <summary>
        /// The coordinate of the point along the single axis Xi. 
        /// </summary>
        public double Xi { get { return xi; } }

        /// <summary>
        /// Vector with the coordinates of the point. Length = 1.
        /// </summary>
        public double[] Coordinates { get { return new double[] { xi }; } }

        public override string ToString()
        {
            return $"(xi) = ({xi})";
        }
    }
}
