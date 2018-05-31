using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Exceptions;
using ISAAR.MSolve.XFEM.Geometry.Mesh;

//TODO: All this logic is common to any FEM like method. Decouple from XFEM (use Vertex and Cell).
namespace ISAAR.MSolve.XFEM.Entities.Decomposition
{
    /// <summary>
    /// Automatic domain decomposition given initial regions as guides. Since the boundaries between these regions may intersect
    /// elements, the final decomposition will have different boundaries.
    /// </summary>
    class GuideDecomposer: IDomainDecomposer
    {
        private readonly int numRegions;
        private readonly IRegion2D[] guides;
        private readonly BiMesh2D mesh;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="guides">They must be convex.</param>
        public GuideDecomposer(IRegion2D[] guides, BiMesh2D mesh)
        {
            this.guides = guides;
            this.numRegions = guides.Length;
            this.mesh = mesh;
        }

        public XCluster2D CreateSubdomains()
        {
            HashSet<XNode2D>[] internalNodes, boundaryNodes;
            (internalNodes, boundaryNodes) = PartitionNodes();
            HashSet<XContinuumElement2D>[] subdomains = PartitionElements(internalNodes, boundaryNodes);
            (internalNodes, boundaryNodes) = RepartitionNodes(subdomains);

            var cluster = new XCluster2D();
            for (int i = 0; i < numRegions; ++i)
            {
                cluster.AddSubdomain(new XSubdomain2D(i, subdomains[i], internalNodes[i], boundaryNodes[i]));
            }
            return cluster;
        }

        public void UpdateSubdomains(XCluster2D cluster)
        {
            // Do nothing
        }

        private int DecideElementRegion(XContinuumElement2D element, HashSet<XNode2D>[] internalNodes,
            HashSet<XNode2D>[] boundaryNodes)
        {
            foreach (var node in element.Nodes)
            {
                for (int i = 0; i < numRegions; ++i)
                {
                    // The first region that contains part of the element is chosen.
                    //TODO: a more sophisticated criterion, such as looking at the relative areas is needed.
                    if (internalNodes[i].Contains(node)) return i;
                }
            }
            throw new IncorrectDecompositionException("This element's nodes do not belong to any of the regions provided");
        }

        private (HashSet<XNode2D>[] internalNodes, HashSet<XNode2D>[] boundaryNodes) PartitionNodes()
        {
            var internalNodes = new HashSet<XNode2D>[numRegions];
            var boundaryNodes = new HashSet<XNode2D>[numRegions];
            for (int i = 0; i < numRegions; ++i)
            {
                internalNodes[i] = new HashSet<XNode2D>();
                boundaryNodes[i] = new HashSet<XNode2D>();
            }

            foreach (var node in mesh.Vertices)
            {
                #region debug
                //if (node.ID == 1250)
                //    Console.WriteLine();
                #endregion

                int multiplicity = 0;
                for (int i = 0; i < numRegions; ++i)
                {
                    NodePosition pos = guides[i].FindRelativePosition(node);
                    if (pos == NodePosition.Internal)
                    {
                        //TODO: this is not foolproof for checking ovelaps. I need to check after going through all regions & nodes
                        if (multiplicity > 0) throw new ArgumentException(
                           $"The regions overlap, as {node} is internal to region {i} but also contained in an earlier one.");
                        
                        internalNodes[i].Add(node);
                        ++multiplicity;
                    }
                    else if (pos == NodePosition.Boundary)
                    {
                        boundaryNodes[i].Add(node);
                        ++multiplicity;
                    }
                }
                if (multiplicity == 0) throw new IncorrectDecompositionException(
                    $"{node} did not belong to any of the regions provided");
            }
            return (internalNodes, boundaryNodes);
        }

        private HashSet<XContinuumElement2D>[] PartitionElements(HashSet<XNode2D>[] internalNodes, 
            HashSet<XNode2D>[] boundaryNodes)
        {
            var subdomains = new HashSet<XContinuumElement2D>[numRegions];
            for (int i = 0; i < numRegions; ++i) subdomains[i] = new HashSet<XContinuumElement2D>();

            foreach (var element in mesh.Cells)
            {
                int elementRegion = DecideElementRegion(element, internalNodes, boundaryNodes);
                subdomains[elementRegion].Add(element);                
            }
            return subdomains;
        }

        private (HashSet<XNode2D>[] internalNodes, HashSet<XNode2D>[] boundaryNodes) RepartitionNodes(
            HashSet<XContinuumElement2D>[] subdomains)
        {
            var internalNodes = new HashSet<XNode2D>[numRegions];
            var boundaryNodes = new HashSet<XNode2D>[numRegions];
            for (int i = 0; i < numRegions; ++i)
            {
                internalNodes[i] = new HashSet<XNode2D>();
                boundaryNodes[i] = new HashSet<XNode2D>();
            }

            foreach (var node in mesh.Vertices)
            {
                // Find out which subdomains it elements belong to
                var nodeSubdomains = new HashSet<int>();
                foreach (var element in mesh.FindElementsWithNode(node))
                {
                    for (int i = 0; i < numRegions; ++i)
                    {
                        if (subdomains[i].Contains(element))
                        {
                            nodeSubdomains.Add(i);
                            break; // Each element only belongs to 1 subdomain exactly.
                        }
                    }
                }

                // Add it to the appropriate subdomans
                if (nodeSubdomains.Count == 1)
                {
                    foreach (var subdomain in nodeSubdomains) internalNodes[subdomain].Add(node);
                }
                else if (nodeSubdomains.Count > 0)
                {
                    foreach (var subdomain in nodeSubdomains) boundaryNodes[subdomain].Add(node);
                }
                else throw new IncorrectDecompositionException($"{node} doesn't belong to any of the identified subdomains");
            }

            return (internalNodes, boundaryNodes);
        }
    }
}
