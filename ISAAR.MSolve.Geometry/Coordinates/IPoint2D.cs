using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Geometry.Coordinates
{
    /// <summary>
    /// Point in a 2-dimensional space. Immutable.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IPoint2D
    {
        /// <summary>
        /// Vector with the coordinates of the point. Length = 2.
        /// </summary>
        double[] Coordinates { get; }

        /// <summary>
        /// The first coordinate of the point.
        /// </summary>
        double X1 { get; }

        /// <summary>
        /// The second coordinate of the point.
        /// </summary>
        double X2 { get; }
    }
}
