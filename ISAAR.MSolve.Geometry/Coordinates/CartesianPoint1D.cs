using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Geometry.Coordinates
{
    /// <summary>
    /// Point in a 1-dimensional cartesian coordinate system. Immutable.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CartesianPoint1D: IPoint1D
    {
        protected readonly double x;

        /// <summary>
        /// Instantiates a <see cref="CartesianPoint1D"/>.
        /// </summary>
        /// <param name="x">The coordinate of the point along the single axis X.</param>
        public CartesianPoint1D(double x)
        {
            this.x = x;
        }

        /// <summary>
        /// The coordinate of the point along the single axis X.
        /// </summary>
        public double X { get { return x; } }

        /// <summary>
        /// Vector with the coordinates of the point. Length = 1.
        /// </summary>
        public double[] Coordinates { get { return new double[] { x }; } }

        public override string ToString()
        {
            return $"(x) = ({x})";
        }
    }
}
