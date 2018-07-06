using ISAAR.MSolve.Geometry.Coordinates;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Discretization.Integration.Points
{
    /// <summary>
    /// Integration point (coordinates & weight) defined in the 1D natural coordinate system of a finite element. Immutable.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class GaussPoint1D: NaturalPoint1D
    {
        /// <summary>
        /// Instantiates a <see cref="GaussPoint1D"/>.
        /// </summary>
        /// <param name="xi">The coordinate of the point along the single axis Xi.</param>
        /// <param name="weight">The weight factor of this integration point.</param>
        public GaussPoint1D(double xi, double weight): base(xi)
        {
            this.Weight = weight;
        }

        /// <summary>
        /// The weight factor of this integration point.
        /// </summary>
        public double Weight { get; }

        public override string ToString()
        {
            return $"xi = {xi} - weight = {Weight}";
        }
    }
}
