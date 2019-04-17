using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Geometry.Coordinates
{
    /// <summary>
    /// Point in a 3-dimensional cartesian coordinate system. It can also represent points in 1-dimensional or 2-dimension 
    /// spaces, but not all coordinates will be used. Immutable
    /// Authors: Serafeim Bakalakos, Dimitris Tsapetis
    /// </summary>
    public class CartesianPoint : IPoint
	{
		protected readonly double x;
		protected readonly double y;
		protected readonly double z;

		/// <summary>
		/// Instantiates a <see cref="CartesianPoint"/>
		/// </summary>
		/// <param name="x">The coordinate of the point along axis X.</param>
		/// <param name="y">The coordinate of the point along axis Y.</param>
		/// <param name="z">The coordinate of the point along axis Z.</param>
		public CartesianPoint(double x, double y = 0.0, double z = 0.0)
		{
			this.x = x;
			this.y = y;
			this.z = z;
		}

		/// <summary>
		/// Instantiates a <see cref="CartesianPoint"/>
		/// </summary>
		/// <param name="coordinates">Vector with the coordinates of the point. Length = 3.</param>
		public CartesianPoint(double[] coordinates)
		{
			this.x = coordinates[0];
			this.y = coordinates[1];
			this.z = coordinates[2];
        }

        /// <summary>
        /// Vector with the coordinates of the point. Length = 3.
        /// </summary>
        public double[] Coordinates => new double[] { x, y, z };

        public double X1 => x;

        public double X2 => y;

        public double X3 => z;

        /// <summary>
        /// The coordinate of the point along axis X.
        /// </summary>
        public double X => x;

		/// <summary>
		/// The coordinate of the point along axis Y.
		/// </summary>
		public double Y => y;

		/// <summary>
		/// The coordinate of the point along axis Z.
		/// </summary>
		public double Z => z;

        /// <summary>
        /// Calculates the Euclidian distance between a <see cref="CartesianPoint"/> named <paramref name="other"/> and this 
        /// one. It will be non negative.
        /// </summary>
        /// <param name="other">The other <see cref="CartesianPoint"/>.</param>
        public double CalculateDistanceFrom(CartesianPoint other) //TODO: this should be implemented for IPoint
        {
            double dx = this.x - other.x;
            double dy = this.y - other.y;
            double dz = this.z - other.z;
            return Math.Sqrt(dx * dx + dy * dy + dz * dz);
        }

        public override string ToString() => $"(x, y, z) = ({x}, {y}, {z})";
    }
}
