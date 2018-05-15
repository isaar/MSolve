using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.LevelSetUpdating
{
    class LevelSetUpdaterStolarska : ILevelSetUpdater
    {
        public LevelSetUpdaterStolarska() { }

        public HashSet<XNode2D> Update(ICartesianPoint2D oldTip, double localGrowthAngle, double growthLength, double dx, double dy, 
            IReadOnlyList<XNode2D> allNodes, ISet<XNode2D> crackBodyNodesAll,
            Dictionary<XNode2D, double> levelSetsBody, Dictionary<XNode2D, double> levelSetsTip)
        {
            double unitDx = dx / growthLength;
            double unitDy = dy / growthLength;
            var newTip = new CartesianPoint2D(oldTip.X + dx, oldTip.Y + dy);
            var newSegment = new DirectedSegment2D(oldTip, newTip);

            var crackBodyNodesModified = new HashSet<XNode2D>();
            foreach (XNode2D node in allNodes)
            {
                // Rotate the ALL tip level sets towards the new tip and then advance them
                double rotatedTipLevelSet = (node.X - oldTip.X) * unitDx + (node.Y - oldTip.Y) * unitDy;
                levelSetsTip[node] = rotatedTipLevelSet - growthLength;

                // Only some body level sets are updated (See Stolarska 2001) 
                // WARNING: this is quite dangerous. If the crack turns/kinks sufficiently over its whole history, then nodes far
                // alway from the tip that should remain unchanged, will get different level sets. A new formula is needed to 
                // avoid messing them up. In a narrow band approach, it may be possible to limit the level set update in nodes
                // near to the old tip the new tip or the inbetween segment.
                if (rotatedTipLevelSet > 0.0) 
                {
                    levelSetsBody[node] = newSegment.SignedDistanceOf(node);
                    if (crackBodyNodesAll.Contains(node)) crackBodyNodesModified.Add(node);
                }
            }

            if (crackBodyNodesModified.Count != 0)
            {
                Console.WriteLine("There are nodes that were previously enriched with Heaviside, but now their level set"
                    + " changes. For 1st order LSM, this usually indicates that the crack has turned sufficiently towards"
                    + " its original geometry. In that case, the body level set of these nodes is computed from the crack"
                    + " extension far away from their position and is thus incorrect. These nodes are: ");
                foreach (var node in crackBodyNodesModified) Console.WriteLine(node);
                Console.WriteLine();
            }
            return crackBodyNodesModified;
        }
    }
}
