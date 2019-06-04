using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Solvers
{
    public class CrackedFetiDPSubdomainCornerNodes : ICornerNodeSelection
    {
        private readonly ICrackDescription crack;
        private readonly Dictionary<int, HashSet<INode>> currentCornerNodes;

        public CrackedFetiDPSubdomainCornerNodes(ICrackDescription crack,
           Dictionary<int, HashSet<INode>> initialCornerNodes)
        {
            this.crack = crack;
            this.currentCornerNodes = initialCornerNodes;
        }

        public Dictionary<int, HashSet<INode>> SelectCornerNodesOfSubdomains()
        {
            HashSet<XNode> enrichedBoundaryNodes = FindNewEnrichedBoundaryNodes();

            foreach (XNode node in enrichedBoundaryNodes)
            {
                foreach (int subdomainID in node.SubdomainsDictionary.Keys)
                {
                    currentCornerNodes[subdomainID].Add(node);
                }
            }
            return currentCornerNodes;
        }

        private HashSet<XNode> FindNewEnrichedBoundaryNodes()
        {
            var enrichedBoundaryNodes = new HashSet<XNode>();
            foreach (CartesianPoint crackTip in crack.CrackTips)
            {
                foreach (XContinuumElement2D tipElement in crack.CrackTipElements[crackTip])
                {
                    foreach (XNode node in tipElement.Nodes)
                    {
                        if (node.SubdomainsDictionary.Count > 1) enrichedBoundaryNodes.Add(node);
                    }
                }
            }
            foreach (var crackBodyNewNodes in crack.CrackBodyNodesNew)
            {
                foreach (XNode node in crackBodyNewNodes.Value)
                {
                    if (node.SubdomainsDictionary.Count > 1) enrichedBoundaryNodes.Add(node);
                }
            }
            return enrichedBoundaryNodes;
        }
    }
}
