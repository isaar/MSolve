using ISAAR.MSolve.Geometry.Coordinates;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Interpolation.Inverse
{
    /// <summary>
    /// Inverse mapping of an isoparametric interpolation, namely from global cartesian to natural (element local) coordinate 
    /// system. In general these are computationally inefficient and inexact, beacuse they use iterative algorithms, with 
    /// the exception of linear interpolations, for which there exist analytic formulas.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IInverseInterpolation2D
    {
        /// <summary>
        /// Finds the coordinates in the natural (element local) system of a point known in the global cartesian system.
        /// </summary>
        /// <param name="point">The coordinates of the point in the global cartesian system. The point must be internal to the 
        ///     element for the method to work, unless otherwise stated by the implementing class.</param>
        /// <returns></returns>
        NaturalPoint2D TransformPointCartesianToNatural(CartesianPoint2D point);
    }
}
