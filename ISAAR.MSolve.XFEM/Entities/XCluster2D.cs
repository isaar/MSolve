using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Entities
{
    /// <summary>
    /// Inter-subdomain communication 
    /// </summary>
    class XCluster2D
    {
        private readonly List<XSubdomain2D> subdomains;

        public XCluster2D()
        {
            this.subdomains = new List<XSubdomain2D>();
        }

        public XClusterDofOrderer DofOrderer { get; private set; }
        public List<XSubdomain2D> Subdomains { get { return subdomains; } }

        public void AddSubdomain(XSubdomain2D subdomain)
        {
            //TODO: all entities should have a builder that uses sets, etc to make sure other entities are unique, etc and the 
            //main collection that will be used by solvers, which must use sorted and efficient collections.
            if (subdomains.Contains(subdomain)) throw new ArgumentException(
                $"A subdomain with id = {subdomain.ID} already exists.");  //TODO: Write a DuplicateEntryException for entities.
            subdomains.Add(subdomain);
        }

        public void AddSubdomains(IEnumerable<XSubdomain2D> subdomains)
        {
            foreach (var subdomain in subdomains)
            {
                AddSubdomain(subdomain);
            }
        }

        //TODO: Should this be a dedicated class?
        //TODO: Should this be done, when a subdomain is added? What if a subdomain is added first and then nodes are added
        //      to the subdomain.
        public Dictionary<XNode2D, SortedSet<XSubdomain2D>> FindBoundaryNodeMembership()
        {
            // Assign subdomains to boundary nodes
            var membership = new Dictionary<XNode2D, SortedSet<XSubdomain2D>>();
            foreach (var subdomain in subdomains)
            {
                foreach (var node in subdomain.BoundaryNodes)
                {
                    bool isTracked = membership.TryGetValue(node, out SortedSet<XSubdomain2D> nodeSubdomains);
                    if (isTracked) nodeSubdomains.Add(subdomain);
                    else
                    {
                        nodeSubdomains = new SortedSet<XSubdomain2D>();
                        nodeSubdomains.Add(subdomain);
                        membership.Add(node, nodeSubdomains);
                    }
                }
            }

            // Check that the multiplicity of each boundary node is at least 2. Otherwise it would have been an internal node.
            foreach (var nodeSubdomains in membership)
            {
                //TODO: dedicated exception class for BadModelConstruction
                if (nodeSubdomains.Value.Count < 2) throw new Exception(
                    $"Found boundary {nodeSubdomains.Key}, with multiplicity = {nodeSubdomains.Value.Count} < 2");
            }
            return membership;
        }

        public Dictionary<XNode2D, SortedSet<XSubdomain2D>> FindEnrichedBoundaryNodeMembership()
        {
            // Assign subdomains to boundary nodes
            var membership = new Dictionary<XNode2D, SortedSet<XSubdomain2D>>();
            foreach (var subdomain in subdomains)
            {
                foreach (var node in subdomain.BoundaryNodes)
                {
                    if (node.EnrichmentItems.Count == 0) continue;
                    bool isTracked = membership.TryGetValue(node, out SortedSet<XSubdomain2D> nodeSubdomains);
                    if (isTracked) nodeSubdomains.Add(subdomain);
                    else
                    {
                        nodeSubdomains = new SortedSet<XSubdomain2D>();
                        nodeSubdomains.Add(subdomain);
                        membership.Add(node, nodeSubdomains);
                    }
                }
            }

            // Check that the multiplicity of each boundary node is at least 2. Otherwise it would have been an internal node.
            foreach (var nodeSubdomains in membership)
            {
                //TODO: dedicated exception class for BadModelConstruction
                if (nodeSubdomains.Value.Count < 2) throw new Exception(
                    $"Found boundary {nodeSubdomains.Key}, with multiplicity = {nodeSubdomains.Value.Count} < 2");
            }
            return membership;
        }

        public SortedSet<XSubdomain2D> FindEnrichedSubdomains()
        {
            var enrichedSubdomains = new SortedSet<XSubdomain2D>();
            foreach (var subdomain in subdomains)
            {
                if (subdomain.HasEnrichedNodes()) enrichedSubdomains.Add(subdomain);
            }
            return enrichedSubdomains;
        }

        public XSubdomain2D FindSubdomainOfElement(XContinuumElement2D element)
        {
            foreach (var subdomain in subdomains)
            {
                if (subdomain.Elements.Contains(element)) return subdomain;
            }
            throw new KeyNotFoundException("This element does not belong to any subdomain.");
        }

        public void OrderStandardDofs(Model2D model)
        {
            DofOrderer = XClusterDofOrderer.CreateNodeMajor(model, this);
        }
    }
}
