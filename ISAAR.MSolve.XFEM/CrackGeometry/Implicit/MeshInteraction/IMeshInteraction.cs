using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Elements;

namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.MeshInteraction
{
    interface IMeshInteraction
    {
        CrackElementPosition FindRelativePositionOf(XContinuumElement2D element);
    }
}
