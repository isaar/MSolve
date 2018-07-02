using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.Discretization.Integration.Points
{
    /// <summary>
    /// Integration point (coordinates & weight) defined in the 2D natural coordinate system of a finite element. Immutable.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class GaussPoint2D : NaturalPoint2D
    {
        /// <summary>
        /// Instantiates a <see cref="GaussPoint2D"/>.
        /// </summary>
        /// <param name="xi">The coordinate of the point along natural axis Xi.</param>
        /// <param name="eta">The coordinate of the point along natural axis Eta.</param>
        /// <param name="weight">The weight factor of this integration point.</param>
        public GaussPoint2D(double xi, double eta, double weight): base (xi, eta)
        {
            this.Weight = weight;
        }

        /// <summary>
        /// The weight factor of this integration point.
        /// </summary>
        public double Weight { get; }

        public override string ToString()
        {
            return $"(xi, eta) = ({xi}, {eta}) - weight = {Weight}";
        }
    }
}
