using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

//TODO: different methods with local polar input or with global cartesian input
namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    class LevelSetUpdaterOLD: ILevelSetUpdater
    {
        public LevelSetUpdaterOLD()
        {
        }

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
                double rotatedTipLevelSet = (node.X - newTip.X) * unitDx + (node.Y - newTip.Y) * unitDy;
                levelSetsTip[node] = rotatedTipLevelSet - newSegment.Length;

                if (rotatedTipLevelSet > 0.0) // Only some body level sets are updated (See Stolarska 2001) 
                {
                    levelSetsBody[node] = newSegment.SignedDistanceOf(node);
                }
            }
        }
    }
}
