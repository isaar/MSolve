using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Exceptions;

namespace ISAAR.MSolve.XFEM.Entities.Decomposition
{
    // Can be abstracted for any vertices, cells
    class NonOverlappingRegionDecomposer2D: IDomainDecomposer
    {
        private readonly IReadOnlyList<XNode2D> nodes;
        private readonly IReadOnlyList<XContinuumElement2D> elements;
        private readonly IReadOnlyList<IRegion2D> regions;

        public NonOverlappingRegionDecomposer2D(IReadOnlyList<XNode2D> nodes, IReadOnlyList<XContinuumElement2D> elements,
            IReadOnlyList<IRegion2D> regions)
        {
            this.nodes = nodes;
            this.elements = elements;
            this.regions = regions;
        }

        public XCluster2D CreateSubdomains()
        {
            XSubdomain2D[] subdomains = AssignNodesToSubdomains(nodes, regions);
            AddElementsToSubdomains(elements, subdomains);
            var cluster = new XCluster2D();
            cluster.AddSubdomains(subdomains);
            return cluster;
        }

        public void UpdateSubdomains(XCluster2D cluster)
        {
            // Do nothing
        }

        private void AddElementsToSubdomains(IReadOnlyList<XContinuumElement2D> elements, XSubdomain2D[] subdomains)
        {
            foreach (var element in elements)
            {
                bool added = false;
                foreach (var subdomain in subdomains)
                {
                    added = subdomain.AddElementIfInternal(element);
                    if (added) break;
                }
                if (!added) throw new IncorrectDecompositionException(
                    "An element did not belong to any of the regions provided");
            }
        }

        private XSubdomain2D[] AssignNodesToSubdomains(IReadOnlyList<XNode2D> nodes, IReadOnlyList<IRegion2D> regions)
        {
            var subdomains = new XSubdomain2D[regions.Count];
            for (int s = 0; s < regions.Count; ++s) subdomains[s] = new XSubdomain2D(s);

            foreach (var node in nodes)
            {
                int multiplicity = 0;
                for (int s = 0; s < regions.Count; ++s)
                {
                    NodePosition pos = regions[s].FindRelativePosition(node);
                    if (pos == NodePosition.Internal)
                    {
                        if (multiplicity > 0) throw new ArgumentException(
                            $"The regions overlap, as {node} is internal to region {s} and boundary to an earlier one.");
                        subdomains[s].AddInternalNode(node);
                        ++multiplicity;
                        break; // It cannot belong to other regions. TODO: perhaps this should be checked
                    }
                    else if (pos == NodePosition.Boundary)
                    {
                        ++multiplicity;
                        subdomains[s].AddBoundaryNode(node);
                        // It can also belong to other regions.
                    }
                }
                if (multiplicity == 0) throw new IncorrectDecompositionException(
                    $"{node} did not belong to any of the regions provided");
            }

            return subdomains;
        }
    }
}
