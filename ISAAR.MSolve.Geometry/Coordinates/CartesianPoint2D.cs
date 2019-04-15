using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

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
        /// Instantiates a <see cref="CartesianPoint2D"/>.
        /// </summary>
        /// <param name="coordinates">Vector with the coordinates of the point. Length = 2.</param>
        public CartesianPoint2D(Vector2 coordinates)
        {
            this.x = coordinates[0];
            this.y = coordinates[1];
        }

        /// <summary>
        /// Vector with the coordinates of the point. Length = 2.
        /// </summary>
        public double[] Coordinates => new double[] { x, y };

        /// <summary>
        /// The coordinate of the point along axis X.
        /// </summary>
        public double X => x;

        /// <summary>
        /// The coordinate of the point along axis Y.
        /// </summary>
        public double Y => y;

        public double X1 => x;

        public double X2 => y;

        /// <summary>
        /// Calculates the Euclidian distance between a <see cref="CartesianPoint2D"/> named <paramref name="other"/> and this 
        /// one. It will be non negative.
        /// </summary>
        /// <param name="other">The other <see cref="CartesianPoint2D"/>.</param>
        public double CalculateDistanceFrom(CartesianPoint2D other) //TODO: this should be implemented for IPoint2D
        {
            double dx = this.x - other.x;
            double dy = this.y - other.y;
            return Math.Sqrt(dx * dx + dy * dy);
        }

        public override string ToString() => $"(x, y) = ({x}, {y})";
    }
}
