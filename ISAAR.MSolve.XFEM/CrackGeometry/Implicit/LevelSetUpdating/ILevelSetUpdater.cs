using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.LevelSetUpdating
{
    interface ILevelSetUpdater
    {
        void Update(ICartesianPoint2D oldTip, double localGrowthAngle, double growthLength, double dx, double dy,
            IReadOnlyList<XNode2D> allNodes, Dictionary<XNode2D, double> levelSetsBody, Dictionary<XNode2D, double> levelSetsTip);
    }
}
