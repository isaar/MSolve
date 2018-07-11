using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Geometry.Coordinates
{
	/// <summary>
	/// Point in a 3-dimensional cartesian coordinate system. Immutable
	/// Authors: Dimitris Tsapetis
	/// </summary>
    public class CartesianPoint3D
	{
		protected readonly double x;
		protected readonly double y;
		protected readonly double z;

		/// <summary>
		/// Instantiates a <see cref="CartesianPoint3D"/>
		/// </summary>
		/// <param name="x">The coordinate of the point along axis X.</param>
		/// <param name="y">The coordinate of the point along axis Y.</param>
		/// <param name="z">The coordinate of the point along axis Z.</param>
		public CartesianPoint3D(double x, double y, double z)
		{
			this.x = x;
			this.y = y;
			this.z = z;
		}

		/// <summary>
		/// Instantiates a <see cref="CartesianPoint3D"/>
		/// </summary>
		/// <param name="coordinates">Vector with the coordinates of the point. Length = 3.</param>
		public CartesianPoint3D(double[] coordinates)
		{
			this.x = coordinates[0];
			this.y = coordinates[1];
			this.z = coordinates[2];
		}

		/// <summary>
		/// The coordinate of the point along axis X.
		/// </summary>
		public double X=>x;

		/// <summary>
		/// The coordinate of the point along axis Y.
		/// </summary>
		public double Y => y;

		/// <summary>
		/// The coordinate of the point along axis Z.
		/// </summary>
		public double Z => z;

		/// <summary>
		/// Vector with the coordinates of the point. Length = 3.
		/// </summary>
		public double[] Coordinates => new double[] {x, y, z};

		public override string ToString()
		{
			return $"(x, y, z) = ({x}, {y}, {z})";
		}
	}
}
