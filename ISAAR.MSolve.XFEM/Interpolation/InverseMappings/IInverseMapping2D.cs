using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.XFEM.Interpolation.InverseMappings
{
    interface IInverseMapping2D
    {
        NaturalPoint2D TransformCartesianToNatural(CartesianPoint2D point);
    }
}
