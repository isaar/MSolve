using System;
using System.Collections.Generic;
using System.Text;
using System.Xml;

namespace ISAAR.MSolve.Geometry.Coordinates
{
	/// <summary>
	/// Point in a 3-dimensional cartesian coordinate system, which is local to the finite element. Immutable
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class NaturalPoint3D : IPoint3D
	{
		protected readonly double xi;
		protected readonly double eta;
		protected readonly double zeta;

		/// <summary>
		/// Instantiates a <see cref="NaturalPoint3D"/>
		/// </summary>
		/// <param name="xi">The coordinate of the point along local axis Xi</param>
		/// <param name="eta">The coordinate of the point along local axis Eta</param>
		/// <param name="zeta">The coordinate of the point along local axis Zeta</param>
		public NaturalPoint3D(double xi, double eta, double zeta)
		{
			this.xi = xi;
			this.eta = eta;
			this.zeta = zeta;
		}

		/// <summary>
		/// Instantiates a <see cref="NaturalPoint3D"/>
		/// </summary>
		/// <param name="coordinates">Vector with the coordinates of the point. Length = 2.</param>
		public NaturalPoint3D(double[] coordinates)
		{
			this.xi = coordinates[0];
			this.eta = coordinates[1];
			this.zeta = coordinates[2];
		}

		/// <summary>
		/// The coordinate of the point along local axis Xi.
		/// </summary>
		public double Xi => xi;

		/// <summary>
		/// The coordinate of the point along local axis Eta.
		/// </summary>
		public double Eta => eta;

		/// <summary>
		/// The coordinate of the point along local axis Zeta.
		/// </summary>
		public double Zeta => zeta;

		/// <summary>
		/// Vector with the coordinates of the point. Length = 3.
		/// </summary>
		public double[] Coordinates => new double[] {xi, eta, zeta};

		public override string ToString()
		{
			return $"(xi, eta, zeta)=({xi}, {eta}, {zeta})";
		}
	}
}
