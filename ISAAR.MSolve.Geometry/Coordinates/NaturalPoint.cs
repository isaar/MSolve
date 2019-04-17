using System;
using System.Collections.Generic;
using System.Text;
using System.Xml;

namespace ISAAR.MSolve.Geometry.Coordinates
{
    /// <summary>
    /// Point in a 3-dimensional cartesian coordinate system, which is local to the finite element. It can also represent points 
    /// in 1-dimensional or 2-dimension spaces, but not all coordinates will be used. Immutable
    /// Authors: Serafeim Bakalakos, Dimitris Tsapetis
    /// </summary>
    public class NaturalPoint : IPoint
	{
		protected readonly double xi;
		protected readonly double eta;
		protected readonly double zeta;

		/// <summary>
		/// Instantiates a <see cref="NaturalPoint"/>
		/// </summary>
		/// <param name="xi">The coordinate of the point along local axis Xi</param>
		/// <param name="eta">The coordinate of the point along local axis Eta</param>
		/// <param name="zeta">The coordinate of the point along local axis Zeta</param>
		public NaturalPoint(double xi, double eta = 0.0, double zeta = 0.0)
		{
			this.xi = xi;
			this.eta = eta;
			this.zeta = zeta;
		}

		/// <summary>
		/// Instantiates a <see cref="NaturalPoint"/>
		/// </summary>
		/// <param name="coordinates">Vector with the coordinates of the point. Length = 2.</param>
		public NaturalPoint(double[] coordinates)
		{
			this.xi = coordinates[0];
			this.eta = coordinates[1];
			this.zeta = coordinates[2];
		}

        public double X1 => xi;

        public double X2 => eta;

        public double X3 => zeta;

        /// <summary>
        /// Vector with the coordinates of the point. Length = 3.
        /// </summary>
        public double[] Coordinates => new double[] { xi, eta, zeta };

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

        public override string ToString() => $"(xi, eta, zeta)=({xi}, {eta}, {zeta})";
    }
}
