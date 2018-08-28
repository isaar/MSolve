using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Geometry.Coordinates
{
    /// <summary>
    /// Point in a 2-dimensional cartesian coordinate system. Immutable.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CartesianPoint2D: IPoint2D
    {
        protected readonly double x;
        protected readonly double y;

        /// <summary>
        /// Instantiates a <see cref="CartesianPoint2D"/>.
        /// </summary>
        /// <param name="x">The coordinate of the point along axis X.</param>
        /// <param name="y">The coordinate of the point along axis Y.</param>
        public CartesianPoint2D(double x, double y)
        {
            this.x = x;
            this.y = y;
        }

        /// <summary>
        /// Instantiates a <see cref="CartesianPoint2D"/>.
        /// </summary>
        /// <param name="coordinates">Vector with the coordinates of the point. Length = 2.</param>
        public CartesianPoint2D(double[] coordinates)
        {
            this.x = coordinates[0];
            this.y = coordinates[1];
        }

        /// <summary>
        /// The coordinate of the point along axis X.
        /// </summary>
        public double X { get { return x; } }

        /// <summary>
        /// The coordinate of the point along axis Y.
        /// </summary>
        public double Y { get { return y; } }

        /// <summary>
        /// Vector with the coordinates of the point. Length = 2.
        /// </summary>
        public double[] Coordinates { get { return new double[] { x, y }; } }

        public override string ToString()
        {
            return $"(x, y) = ({x}, {y})";
        }
    }
}
