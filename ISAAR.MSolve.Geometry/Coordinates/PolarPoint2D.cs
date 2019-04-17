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
    public class PolarPoint2D : IPoint
    {
        protected readonly double r;
        protected readonly double theta;

        /// <summary>
        /// Instantiates a <see cref="PolarPoint2D"/>.
        /// </summary>
        /// <param name="r">The radial coordinate of the point.</param>
        /// <param name="theta">The angular coordinate of the point.</param>
        public PolarPoint2D(double r, double theta)
        {
            this.r = r;
            this.theta = theta;
        }

        /// <summary>
        /// Vector with the coordinates of the point. Length = 2.
        /// </summary>
        public double[] Coordinates => new double[] { r, theta };

        /// <summary>
        /// The radial coordinate of the point.
        /// </summary>
        public double R => r;

        /// <summary>
        /// The angular coordinate of the point.
        /// </summary>
        public double Theta => theta;

        public double X1 => r;

        public double X2 => theta;

        public double X3 => throw new NotImplementedException();

        public override string ToString() => $"(r, theta) = ({r}, {theta})";
    }
}
