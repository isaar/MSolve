using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.Discretization.Integration
{
    /// <summary>
    /// Integration point (coordinates & weight) defined in the 3D natural coordinate system of a finite element. It can also 
    /// represent points in 1-dimensional or 2-dimension spaces, but not all coordinates will be used. Immutable.
    /// Authors: Serafeim Bakalakos, Dimitris Tsapetis
    /// </summary>
    public class GaussPoint : NaturalPoint
	{
        /// <summary>
        /// Creates a new instance of <see cref="GaussPoint"/> for an 1D coordinate system.
        /// </summary>
        /// <param name="xi"></param>
        /// <param name="weight"></param>
        public GaussPoint(double xi, double weight): base(xi)
        {
            this.Weight = weight;
        }

        /// <summary>
        /// Creates a new instance of <see cref="GaussPoint"/> for a 2D coordinate system.
        /// </summary>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <param name="weight"></param>
        public GaussPoint(double xi, double eta, double weight) : base(xi, eta)
        {
            this.Weight = weight;
        }

        /// <summary>
        /// Creates a new instance of <see cref="GaussPoint"/> for a 3D coordinate system.
        /// </summary>
        /// <param name="xi">The coordinate of the point along natural axis Xi.</param>
        /// <param name="eta">The coordinate of the point along natural axis Eta.</param>
        /// <param name="zeta">The coordinate of the point along natural axis Zeta.</param>
        /// <param name="weight">The weight factor of this integration point.</param>
        public GaussPoint(double xi, double eta, double zeta, double weight) : base(xi, eta, zeta)
		{
			this.Weight = weight;
		}

        /// <summary>
        /// The weight factor of this integration point.
        /// </summary>
		public double Weight { get; }
    }
}
