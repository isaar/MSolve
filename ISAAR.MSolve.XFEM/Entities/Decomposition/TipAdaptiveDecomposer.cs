using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Exceptions;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Entities.Decomposition
{
    class TipAdaptiveDecomposer: IDomainDecomposer //TODO: extend it to more than one tips/cracks
    {
        private readonly ICrackDescription crack;
        private readonly BiMesh2D mesh;
        private readonly IReadOnlyList<IRegion2D> regions;
        private IDomainDecomposer initialDecomposer; //TODO: this should be discarded once used.

        public TipAdaptiveDecomposer(BiMesh2D mesh, IReadOnlyList<IRegion2D> regions, ICrackDescription crack,
            IDomainDecomposer initialDecomposer)
        {
            this.mesh = mesh;
            this.regions = regions;
            this.crack = crack;
            this.initialDecomposer = initialDecomposer;
        }

        public XCluster2D CreateSubdomains()
        {
            XCluster2D cluster = initialDecomposer.CreateSubdomains();
            initialDecomposer = null; // no longer needed and it might take up significant memory for large meshes.
            return cluster;
        }

        public void UpdateSubdomains(XCluster2D cluster)
        {
            //TODO: this can be done without accessing the single crack. We just need to access the same tip element and nodes
            foreach (ISingleCrack singleCrack in crack.SingleCracks)
            {
                // Find the subdomain that contains the crack tip
                ICartesianPoint2D tip = singleCrack.CrackTips[0];
                XContinuumElement2D tipElement = singleCrack.CrackTipElements[tip][0]; // If there are more neighboring ones, we don't need them
                XSubdomain2D tipSubdomain = cluster.FindSubdomainOfElement(tipElement);

                // Find all elements that have at least one tip enriched node
                //TODO: this can be merged with the next subtask
                var tipEnrichedElements = new HashSet<XContinuumElement2D>();
                ISet<XNode2D> tipNodes = singleCrack.CrackTipNodesNew[singleCrack.CrackTipEnrichments];
                foreach (XNode2D node in tipNodes) tipEnrichedElements.UnionWith(mesh.FindElementsWithNode(node));

                // Find which of these elements must be removed and from which subdomains
                var removedElements = new Dictionary<XContinuumElement2D, XSubdomain2D>();
                foreach (var element in tipEnrichedElements)
                {
                    if (!tipSubdomain.Elements.Contains(element))
                    {
                        XSubdomain2D previousSubdomain = cluster.FindSubdomainOfElement(element);
                        removedElements.Add(element, previousSubdomain);
                    }
                }
                if (removedElements.Count == 0) continue;

                // Find the old subdomains of the nodes of the elements to be removed, before moving the elements
                Dictionary<XNode2D, HashSet<XSubdomain2D>> oldNodeMembership = FindOldNodeMembership(cluster, removedElements);

                // Move these elements to the subdomain that contains the crack tip
                //var modifiedSubdomains = new HashSet<XSubdomain2D>();
                foreach (var elementSubdomainPair in removedElements)
                {
                    XContinuumElement2D element = elementSubdomainPair.Key;
                    XSubdomain2D oldSubdomain = elementSubdomainPair.Value;
                    oldSubdomain.Elements.Remove(element);
                    tipSubdomain.Elements.Add(element);
                    //modifiedSubdomains.Add(previousSubdomain);
                    //modifiedSubdomains.Add(tipSubdomain); //TODO: this only needs to be added once
                }

                // Move their nodes to their new subdomains, which also updates the boundaries.
                RemoveNodesFromSubdomains(oldNodeMembership);
                Dictionary<XNode2D, HashSet<XSubdomain2D>> newNodeMembership =
                    FindNewNodeMembership(cluster, oldNodeMembership.Keys);
                AddNodesToSubdomains(newNodeMembership);
            }
        }

        // TODO: perhaps adding and deleting must be done simultaneously with building the node memberships
        private static void AddNodesToSubdomains(Dictionary<XNode2D, HashSet<XSubdomain2D>> newNodeMembership)
        {
            foreach (var nodeSubdomains in newNodeMembership)
            {
                XNode2D node = nodeSubdomains.Key;
                int nodeSubdomainsCount = nodeSubdomains.Value.Count;
                if (nodeSubdomainsCount > 1) // Will be boundary node between 2 or more subdomains
                {
                    foreach (var subdomain in nodeSubdomains.Value) subdomain.AddBoundaryNode(node);
                }
                else if (nodeSubdomainsCount == 1) // Will be internal to 1 subdomain
                {
                    foreach (var subdomain in nodeSubdomains.Value) subdomain.AddInternalNode(node);
                }
                else throw new IncorrectDecompositionException(
                    $"{nodeSubdomains.Key} didn't belong to any of the identified subdomains");
            }
        }

        private static void RemoveNodesFromSubdomains(Dictionary<XNode2D, HashSet<XSubdomain2D>> oldNodeMembership)
        {
            foreach (var nodeSubdomains in oldNodeMembership)
            {
                XNode2D node = nodeSubdomains.Key;
                int nodeSubdomainsCount = nodeSubdomains.Value.Count;
                if (nodeSubdomainsCount > 1) // Was boundary node between 2 or more subdomains
                {
                    foreach (var subdomain in nodeSubdomains.Value)
                    {
                        subdomain.BoundaryNodes.Remove(node);
                        subdomain.AllNodes.Remove(node);
                    }
                }
                else if (nodeSubdomainsCount == 1) // Was internal to 1 subdomain
                {
                    foreach (var subdomain in nodeSubdomains.Value)
                    {
                        subdomain.InternalNodes.Remove(node);
                        subdomain.AllNodes.Remove(node);
                    }
                }
                else throw new IncorrectDecompositionException(
                    $"{nodeSubdomains.Key} didn't belong to any of the identified subdomains");
            }
        }

        //TODO: perhaps the node membership should be stored in Cluster for all nodes of the model and updated whenever sth changes 
        private Dictionary<XNode2D, HashSet<XSubdomain2D>> FindOldNodeMembership(XCluster2D cluster,
            Dictionary<XContinuumElement2D, XSubdomain2D> removedElements)
        {
            var nodeMembership = new Dictionary<XNode2D, HashSet<XSubdomain2D>>();
            foreach (var elementSubdomainPair in removedElements)
            {
                foreach (var node in elementSubdomainPair.Key.Nodes)
                {
                    if (nodeMembership.ContainsKey(node)) continue; // Avoid nodes that were processed in previous elements
                    var nodeSubdomains = new HashSet<XSubdomain2D>();
                    foreach (var nodeElement in mesh.FindElementsWithNode(node))
                    {
                        // Find the subdomain of this element efficiently
                        //TODO: cache the searched elements for other nodes
                        bool isRemoved = removedElements.TryGetValue(nodeElement, out XSubdomain2D previousSubdomain);
                        if (isRemoved) nodeSubdomains.Add(previousSubdomain);
                        else nodeSubdomains.Add(cluster.FindSubdomainOfElement(nodeElement)); 
                    }
                    nodeMembership.Add(node, nodeSubdomains);
                }
            }
            return nodeMembership;
        }

        private Dictionary<XNode2D, HashSet<XSubdomain2D>> FindNewNodeMembership(XCluster2D cluster, IEnumerable<XNode2D> nodes)
        {
            var nodeMembership = new Dictionary<XNode2D, HashSet<XSubdomain2D>>();
            var elementMembership = new Dictionary<XContinuumElement2D, XSubdomain2D>(); // cache each element for all its nodes
            foreach (var node in nodes)
            {
                var nodeSubdomains = new HashSet<XSubdomain2D>();
                foreach (var element in mesh.FindElementsWithNode(node))
                {
                    // Find the subdomain of this element efficiently
                    bool isProccessed = elementMembership.TryGetValue(element, out XSubdomain2D subdomain);
                    if (!isProccessed)
                    {
                        subdomain = cluster.FindSubdomainOfElement(element);
                        elementMembership.Add(element, subdomain);
                    }
                    nodeSubdomains.Add(subdomain);
                }
                nodeMembership.Add(node, nodeSubdomains);
            }
            return nodeMembership;
        }
    }
}
