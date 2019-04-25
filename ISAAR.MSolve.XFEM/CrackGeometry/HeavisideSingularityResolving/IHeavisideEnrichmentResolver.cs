using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.CrackGeometry.HeavisideSingularityResolving
{
    public interface IHeavisideSingularityResolver
    {
        ISet<XNode> FindHeavisideNodesToRemove(ISingleCrack crack, IMesh2D<XNode, XContinuumElement2D> mesh, 
            ISet<XNode> heavisideNodes);

        ISet<XNode> FindHeavisideNodesToRemove(ISingleCrack crack, IReadOnlyList<XNode> heavisideNodes,
            IReadOnlyList<ISet<XContinuumElement2D>> nodalSupports);
    }
}
