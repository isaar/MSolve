using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

//TODO: should it just pull properties out of the LSM rather than all these parameters? It would be much more abstracted.
//TODO: this strategy interface should probably handle all geometry updates, not only the level sets.
namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.LevelSetUpdating
{
    interface ILevelSetUpdater
    {
        /// <summary>
        /// Returns nodes that were previously enriched with Heaviside and now their body level set changes.
        /// </summary>
        /// <returns></returns>
        HashSet<XNode2D> Update(ICartesianPoint2D oldTip, double localGrowthAngle, double growthLength, double dx, double dy,
            IReadOnlyList<XNode2D> allNodes, ISet<XNode2D> crackBodyNodesAll, 
            Dictionary<XNode2D, double> levelSetsBody, Dictionary<XNode2D, double> levelSetsTip);
    }
}
