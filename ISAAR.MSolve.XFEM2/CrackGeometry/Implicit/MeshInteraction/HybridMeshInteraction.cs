using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.MeshInteraction
{
    class HybridMeshInteraction : IMeshInteraction
    {
        private readonly TrackingExteriorCrackLSM lsm;

        public HybridMeshInteraction(TrackingExteriorCrackLSM lsm)
        {
            this.lsm = lsm;
        }

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

            //Warning: this might actually be worse than Stolarska's criterion. At least that one enriched the dubious 
            //intersected elements with tip enrichments. This one just ignores them.
            if (minBodyLevelSet * maxBodyLevelSet <= 0.0)
            {
                if (minTipLevelSet * maxTipLevelSet <= 0)
                {
                    var outline = ConvexPolygon2D.CreateUnsafe(element.Nodes);
                    if (outline.IsPointInsidePolygon(crackTip)) return CrackElementPosition.ContainsTip;
                }
                else if (maxTipLevelSet < 0) return CrackElementPosition.Intersected;
            }
            return CrackElementPosition.Irrelevant;
        }
    }
}
