using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

//TODO: Delete this after testing. The new version is more efficient and does the same.
namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.MeshInteraction
{
    class HybridMeshInteractionOLD: IMeshInteraction
    {
        private readonly TrackingExteriorCrackLSM lsm;

        public HybridMeshInteractionOLD(TrackingExteriorCrackLSM lsm)
        {
            this.lsm = lsm;
        }

        //TODO: replace this with the Faster~() version
        public CrackElementPosition FindRelativePositionOf(XContinuumElement2D element)
        {
            ICartesianPoint2D crackTip = lsm.GetCrackTip(CrackTipPosition.Single);
            double minBodyLevelSet = double.MaxValue;
            double maxBodyLevelSet = double.MinValue;
            double minTipLevelSet = double.MaxValue;
            double maxTipLevelSet = double.MinValue;

            foreach (XNode2D node in element.Nodes)
            {
                double bodyLevelSet = lsm.LevelSetsBody[node];
                double tipLevelSet = lsm.LevelSetsTip[node];
                if (bodyLevelSet < minBodyLevelSet) minBodyLevelSet = bodyLevelSet;
                if (bodyLevelSet > maxBodyLevelSet) maxBodyLevelSet = bodyLevelSet;
                if (tipLevelSet < minTipLevelSet) minTipLevelSet = tipLevelSet;
                if (tipLevelSet > maxTipLevelSet) maxTipLevelSet = tipLevelSet;
            }

            // Warning: This criterion might give false positives for tip elements (see Serafeim's thesis for details)
            if (minBodyLevelSet * maxBodyLevelSet <= 0.0)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Nodes);
                if (outline.IsPointInsidePolygon(crackTip)) return CrackElementPosition.ContainsTip;
                else if (maxTipLevelSet < 0) return CrackElementPosition.Intersected;
            }
            return CrackElementPosition.Irrelevant;
        }
    }
}
