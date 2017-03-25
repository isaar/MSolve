using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Materials.Crack
{
    interface ICrackMaterialFactory2D
    {
        CrackMaterial2D FindMaterialAtPoint(ICartesianPoint2D point);
    }
}
