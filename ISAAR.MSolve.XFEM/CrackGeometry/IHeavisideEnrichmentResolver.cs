using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.Mesh;

namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    interface IHeavisideEnrichmentResolver
    {
        ISet<XNode2D> FindHeavisideNodesToRemove(IMesh2D<XNode2D, XContinuumElement2D> mesh, ISet<XNode2D> heavisideNodes);
    }
}
