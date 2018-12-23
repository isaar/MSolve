using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Interpolation.InverseMappings
{
    interface IInverseMapping2D
    {
        INaturalPoint2D TransformCartesianToNatural(ICartesianPoint2D point);
    }
}
