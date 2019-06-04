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

//TODO: Do not delete this class until investigating theoritacally why exactly it does not work.
//TODO: In the original class that does work, try to only set the enriched dofs as corner dofs, instead of the whole boundary 
//      enriched nodes
namespace ISAAR.MSolve.XFEM.Solvers
{
    /// <summary>
    /// In <see cref="CrackedFetiDPSubdomainCornerNodes"/> I set as corner nodes all boundary enriched nodes, which happens when 
    /// the crack completely intersect a subdomain. Here I tried to improve upon that idea by only setting as corner nodes the
    /// boundary enriched nodes that lie closest to the crack; one in the positive area and one in the negative. Interestingly 
    /// this does not work. It seems that if there are two boundary Heaviside nodes with positive signed distance, they both
    /// introduce rigid body motions and should be set as corner nodes.
    /// </summary>
    public class CrackedFetiDPSubdomainCornerNodes_WRONG //: ICornerNodeSelection
    {
        private readonly ICrackDescription crack;
        private readonly Dictionary<int, HashSet<INode>> currentCornerNodes;

        public CrackedFetiDPSubdomainCornerNodes_WRONG(ICrackDescription crack,
           Dictionary<int, HashSet<INode>> initialCornerNodes)
        {
            this.crack = crack;
            this.currentCornerNodes = initialCornerNodes;
        }

        public Dictionary<int, HashSet<INode>> SelectCornerNodesOfSubdomains()
        {
            HashSet<XNode> enrichedBoundaryNodes = FindBoundaryNodesOnEachSideOfTheCrack();
            foreach (XNode node in enrichedBoundaryNodes)
            {
                foreach (int subdomainID in node.SubdomainsDictionary.Keys)
                {
                    currentCornerNodes[subdomainID].Add(node);
                }
            }
            return currentCornerNodes;
        }

        private HashSet<XNode> FindBoundaryNodesOnEachSideOfTheCrack()
        {
            var allEnrichedBoundaryNodes = new HashSet<XNode>();

            foreach (ISingleCrack singleCrack in crack.SingleCracks)
            {
                var enrichedNodes = new HashSet<XNode>();

                // Find boundary nodes of the tip element
                foreach (CartesianPoint crackTip in singleCrack.CrackTips)
                {
                    foreach (XContinuumElement2D tipElement in singleCrack.CrackTipElements[crackTip])
                    {
                        foreach (XNode node in tipElement.Nodes)
                        {
                            if (node.SubdomainsDictionary.Count > 1) enrichedNodes.Add(node);
                        }
                    }
                }

                // Find boundary nodes of intersected elements. Limit the search to newly enriched nodes.
                foreach (var newHeavisideNodes in singleCrack.CrackBodyNodesNew)
                {
                    foreach (XNode node in newHeavisideNodes.Value)
                    {
                        if (node.SubdomainsDictionary.Count > 1) enrichedNodes.Add(node);
                    }
                }

                // It is possible that there are two boundary enriched nodes in an element where the crack instersects the
                // inter-subdomain boundary. Keep the closer one to the crack in the positive area and the closer one in the 
                // negative area. 
                if (enrichedNodes.Count > 2)
                {
                    // We cannot just check the signed distances of all the nodes gathered so far for this crack. If the crack
                    // grows from 2 tips, then we need to consider the nodes gathered from each tip separately. Therefore we
                    // cluster the nodes we gathered, using their elements.
                    Dictionary<IElement, HashSet<XNode>> commonElements = FindCommonElements(enrichedNodes);
                    foreach (IElement element in commonElements.Keys)
                    {
                        List<XNode> closestEnrichedBoundaryNodes =
                            FindClosestEnrichedBoundaryNodes(singleCrack, element, commonElements[element]);
                        allEnrichedBoundaryNodes.UnionWith(closestEnrichedBoundaryNodes);
                        Debug.Assert(closestEnrichedBoundaryNodes.Count >= 1);
                        Debug.Assert(closestEnrichedBoundaryNodes.Count <= 2);
                    }
                }
                else allEnrichedBoundaryNodes.UnionWith(enrichedNodes);
            }

            return allEnrichedBoundaryNodes;
        }
        private List<XNode> FindClosestEnrichedBoundaryNodes(ISingleCrack crack, IElement element,
            HashSet<XNode> enrichedBoundaryNodes)
        {
            XNode closestPositiveNode = null;
            double minPositiveDistance = double.MaxValue;
            XNode closestNegativeNode = null;
            double maxNegativeDistance = double.MinValue;
            foreach (XNode node in enrichedBoundaryNodes)
            {
                double signedDistance = crack.SignedDistanceOf(node);
                if (signedDistance >= 0)
                {
                    if (signedDistance < minPositiveDistance)
                    {
                        minPositiveDistance = signedDistance;
                        closestPositiveNode = node;
                    }
                }
                else
                {
                    if (signedDistance > maxNegativeDistance)
                    {
                        maxNegativeDistance = signedDistance;
                        closestNegativeNode = node;
                    }
                }
            }

            // It is possible that there is no positive or no negative enriched boundary node for this element. E.g. This could
            //happen if one of them was too far away from the crack and ended up not getting enriched with Heaviside.
            var result = new List<XNode>();
            if (closestPositiveNode != null) result.Add(closestPositiveNode);
            if (closestNegativeNode != null) result.Add(closestNegativeNode);
            return result;
        }

        private Dictionary<IElement, HashSet<XNode>> FindCommonElements(HashSet<XNode> nodes)
        {
            var commonElements = new Dictionary<IElement, HashSet<XNode>>();
            foreach (XNode node in nodes)
            {
                foreach (IElement element in node.ElementsDictionary.Values)
                {
                    foreach (XNode otherNode in element.Nodes)
                    {
                        if (otherNode == node) continue;
                        if (nodes.Contains(otherNode))
                        {
                            bool exists = commonElements.TryGetValue(element, out HashSet<XNode> elementNodes);
                            if (!exists)
                            {
                                elementNodes = new HashSet<XNode>();
                                commonElements[element] = elementNodes;
                            }
                            elementNodes.Add(node);
                            elementNodes.Add(otherNode);
                        }
                    }
                }
            }
            return commonElements;
        }
    }
}
