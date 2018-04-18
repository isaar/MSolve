using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Geometry.Shapes
{
    interface IOpenCurve2D: ICurve2D
    {
        ICartesianPoint2D Start { get; }
        ICartesianPoint2D End { get; }

        /// <summary>
        /// Unit vector. It will coincide with the normal vector if rotated -PI/2.
        /// </summary>
        Tuple<double, double> TangentAtStart { get; }

        /// <summary>
        /// Unit vector. It will coincide with the normal vector if rotated PI/2.
        /// </summary>
        Tuple<double, double> TangentAtEnd { get; }
    }
}
