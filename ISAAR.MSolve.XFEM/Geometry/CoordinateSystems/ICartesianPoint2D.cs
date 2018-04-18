using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//TODO: Replace double[] with Vector
namespace ISAAR.MSolve.XFEM.Geometry.CoordinateSystems
{
    interface ICartesianPoint2D
    {
        double X { get; }
        double Y { get; }

        /// <summary>
        /// Returns an array of length = 2 containing the coordinates. The point is not mutated if the array is changed.
        /// </summary>
        double[] Coordinates { get; }
    }
}
