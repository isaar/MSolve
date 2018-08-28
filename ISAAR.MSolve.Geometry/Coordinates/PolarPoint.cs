using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Geometry.Coordinates
{
    /// <summary>
    /// Point in a 2-dimensional polar coordinate system. Immutable.
    /// For more see https://en.wikipedia.org/wiki/Polar_coordinate_system.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PolarPoint : IPoint2D
    {
        protected readonly double r;
        protected readonly double theta;

        /// <summary>
        /// Instantiates a <see cref="PolarPoint"/>.
        /// </summary>
        /// <param name="r">The radial coordinate of the point.</param>
        /// <param name="theta">The angular coordinate of the point.</param>
        public PolarPoint(double r, double theta)
        {
            this.r = r;
            this.theta = theta;
        }

        /// <summary>
        /// The radial coordinate of the point.
        /// </summary>
        public double R { get { return r; } }

        /// <summary>
        /// The angular coordinate of the point.
        /// </summary>
        public double Theta { get { return theta; } }

        /// <summary>
        /// Vector with the coordinates of the point. Length = 2.
        /// </summary>
        public double[] Coordinates { get { return new double[] { r, theta }; } }

        public override string ToString()
        {
            return $"(r, theta) = ({r}, {theta})";
        }
    }
}
