using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.MeshInteraction
{
    /// <summary>
    /// Correctly identifies the elements around the crack tip that Stolarska's criterion always marks as tip elements. Uses the
    /// intersections of the suspected tip element with the 1st order zero crack body level set. As this is done for only a 
    /// handful elements (or none), the performance hit is negligible. It may also be possible to extend this criterion to 
    /// higher order LSM by using an appropriate crack-element intersection formula.
    /// </summary>
    class SerafeimMeshInteraction : IMeshInteraction
    {
        private readonly TrackingExteriorCrackLSM lsm;

        public SerafeimMeshInteraction(TrackingExteriorCrackLSM lsm)
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

            if (minBodyLevelSet * maxBodyLevelSet > 0.0) return CrackElementPosition.Irrelevant;
            else // The element is intersected by the zero body level set.
            {
                if (minTipLevelSet > 0) return CrackElementPosition.Irrelevant; // intersected by the crack's extension
                else if (maxTipLevelSet < 0) return CrackElementPosition.Intersected;
                else // Stolarska's criterion marks all the next as tip elements
                {
                    Dictionary<ICartesianPoint2D, double> intersections = FindIntersectionsAndTipLevelSets(element);
                    Debug.Assert((intersections.Count == 2) || (intersections.Count == 1)); // 1 is veeeeery improbable
                    double tipLevelSetInter1 = intersections.First().Value;
                    double tipLevelSetInter2 = intersections.Last().Value;
                    if ((tipLevelSetInter1 > 0) && (tipLevelSetInter2 > 0)) return CrackElementPosition.Irrelevant;
                    else if ((tipLevelSetInter1 < 0) && (tipLevelSetInter2 < 0)) return CrackElementPosition.Intersected;
                    else return CrackElementPosition.ContainsTip;
                }
            }
        }

        private Dictionary<ICartesianPoint2D, double> FindIntersectionsAndTipLevelSets(XContinuumElement2D element)
        {
            // Find intersections of element with the zero body level set. 
            //TODO: abstract this procedure and reuse it here and in LSM
            var intersections = new Dictionary<ICartesianPoint2D, double>();
            var nodes = element.Nodes;
            for (int i = 0; i < nodes.Count; ++i)
            {
                XNode2D node1 = nodes[i];
                XNode2D node2 = nodes[(i + 1) % nodes.Count];
                double bodyLevelSet1 = lsm.LevelSetsBody[node1];
                double bodyLevelSet2 = lsm.LevelSetsBody[node2];

                if (bodyLevelSet1 * bodyLevelSet2 < 0.0)
                {
                    // The intersection point between these nodes can be found using the linear interpolation, see 
                    // Sukumar 2001
                    double k = -bodyLevelSet1 / (bodyLevelSet2 - bodyLevelSet1);
                    double x = node1.X + k * (node2.X - node1.X);
                    double y = node1.Y + k * (node2.Y - node1.Y);

                    double tipLevelSet1 = lsm.LevelSetsTip[node1];
                    double tipLevelSet2 = lsm.LevelSetsTip[node2];
                    double tipLevelSet = tipLevelSet1 + k * (tipLevelSet2 - tipLevelSet1);

                    intersections.Add(new CartesianPoint2D(x, y), tipLevelSet);
                }
                else if (bodyLevelSet1 == 0.0) intersections[node1] = lsm.LevelSetsTip[node1]; // TODO: perhaps some tolerance is needed.
                else if (bodyLevelSet2 == 0.0) intersections[node2] = lsm.LevelSetsTip[node2];
            }

            return intersections;
        }
    }
}
