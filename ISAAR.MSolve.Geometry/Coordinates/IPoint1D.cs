using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Geometry.Coordinates
{
    /// <summary>
    /// Point in a 1-dimensional space. Immutable.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IPoint1D
    {
        /// <summary>
        /// Vector with the coordinates of the point. Length = 1.
        /// </summary>
        double[] Coordinates { get; }
    }
}
