using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.Discretization.Integration.Points
{
	/// <summary>
	/// Integration point (coordinates & weight) defined in the 3D natural coordinate system of a finite element. Immutable
	/// Authors: Dimitris Tsapetis
	/// </summary>
    public class GaussPoint3D: NaturalPoint3D
	{
		/// <summary>
		/// Instantiates a <see cref="GaussPoint3D"/>
		/// </summary>
		/// <param name="xi">The coordinate of the point along natural axis Xi.</param>
		/// <param name="eta">The coordinate of the point along natural axis Eta.</param>
		/// <param name="zeta">The coordinate of the point along natural axis Zeta.</param>
		/// <param name="weight">The weight factor of this integration point.</param>
		public GaussPoint3D(double xi, double eta, double zeta, double weight):base(xi,eta,zeta)
		{
			this.Weight = weight;
		}

		public double Weight { get; }
    }
}
