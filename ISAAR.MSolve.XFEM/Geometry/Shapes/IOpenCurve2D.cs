using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.XFEM.Geometry.Shapes
{
    interface IOpenCurve2D: ICurve2D
    {
        CartesianPoint2D Start { get; }
        CartesianPoint2D End { get; }

        /// <summary>
        /// Unit vector. It will coincide with the normal vector if rotated -PI/2.
        /// </summary>
        Vector2 TangentAtStart { get; }

        /// <summary>
        /// Unit vector. It will coincide with the normal vector if rotated PI/2.
        /// </summary>
        Vector2 TangentAtEnd { get; }
    }
}
