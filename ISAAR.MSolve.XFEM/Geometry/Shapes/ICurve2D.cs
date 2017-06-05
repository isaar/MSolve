using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Geometry.Shapes
{
    interface ICurve2D
    {
        double SignedDistance(ICartesianPoint2D point);

        /// <summary>
        /// Unit length
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        Tuple<double, double> NormalThrough(ICartesianPoint2D point);
    }
}
