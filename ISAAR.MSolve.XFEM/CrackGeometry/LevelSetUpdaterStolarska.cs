using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    class LevelSetUpdaterStolarska : ILevelSetUpdater
    {
        public LevelSetUpdaterStolarska() { }

        public void Update(ICartesianPoint2D oldTip, double localGrowthAngle, double growthLength, double dx, double dy, 
            IReadOnlyList<XNode2D> allNodes, Dictionary<XNode2D, double> levelSetsBody, Dictionary<XNode2D, double> levelSetsTip)
        {
            double unitDx = dx / growthLength;
            double unitDy = dy / growthLength;
            var newTip = new CartesianPoint2D(oldTip.X + dx, oldTip.Y + dy);
            var newSegment = new DirectedSegment2D(oldTip, newTip);

            foreach (XNode2D node in allNodes)
            {
                // Rotate the ALL tip level sets towards the new tip and then advance them
                double rotatedTipLevelSet = (node.X - oldTip.X) * unitDx + (node.Y - oldTip.Y) * unitDy;
                levelSetsTip[node] = rotatedTipLevelSet - growthLength;

                if (rotatedTipLevelSet > 0.0) // Only some body level sets are updated (See Stolarska 2001) 
                {
                    levelSetsBody[node] = newSegment.SignedDistanceOf(node);
                }
            }
        }
    }
}
