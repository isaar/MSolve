using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Replace double[] with Vector
namespace ISAAR.MSolve.XFEM.Geometry.CoordinateSystems
{
    interface ICartesianPoint2D
    {
        double X { get; }
        double Y { get; }

        /// <summary>
        /// Returns a vector of length = 2 containing the coordinates. The point is not mutated if the vector is changed.
        /// </summary>
        Vector2 Coordinates { get; }
    }
}
