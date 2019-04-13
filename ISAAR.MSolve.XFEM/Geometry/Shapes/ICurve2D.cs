using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.XFEM.Geometry.Shapes
{
    interface ICurve2D
    {
        double SignedDistance(CartesianPoint2D point);

        /// <summary>
        /// Unit length
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        Vector2 NormalThrough(CartesianPoint2D point);
    }
}
