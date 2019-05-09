using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

//TODO: Use a SortedSuperSet custom set class, to include the internal and boundaries in a superset
//TODO: During construction, I need it to be mutable, but afterwards not. This goes for Model and CLuster as weel. Use views.
//TODO: should these node sets exposed be ISet<>? I depend on their lookup/insertion performance during construction and order
//      during dof ordering.
namespace ISAAR.MSolve.XFEM.Entities
{
    /// <summary>
    /// Only for enriched dofs.
    /// </summary>
    class XSubdomain2D_old: IComparable<XSubdomain2D_old>
    {
        private readonly HashSet<XContinuumElement2D> elements;
        private readonly SortedSet<XNode> allNodes; // Having it sorted is better for ordering
        private readonly HashSet<XNode> boundaryNodes;
        private readonly HashSet<XNode> internalNodes;

        public XSubdomain2D_old(int id)
        {
            this.ID = id;
            this.allNodes = new SortedSet<XNode>();
            this.elements = new HashSet<XContinuumElement2D>();
            this.boundaryNodes = new HashSet<XNode>();
            this.internalNodes = new HashSet<XNode>();
        }

        public XSubdomain2D_old(int id, HashSet<XContinuumElement2D> elements, HashSet<XNode> internalNodes, 
            HashSet<XNode> boundaryNodes)
        {
            this.ID = id;
            this.elements = elements;
            this.internalNodes = internalNodes;
            this.boundaryNodes = boundaryNodes;
            this.allNodes = new SortedSet<XNode>(internalNodes);
            this.allNodes.UnionWith(boundaryNodes);
        }

        public int ID { get; }
        public ISet<XContinuumElement2D> Elements { get { return elements; } }
        public ISet<XNode> AllNodes { get { return allNodes; } }
        public ISet<XNode> BoundaryNodes { get { return boundaryNodes; } } 
        public ISet<XNode> InternalNodes { get { return internalNodes; } }
        public XSubdomainDofOrderer DofOrderer { get; set; }

        public void AddBoundaryNode(XNode node)
        { //TODO: perhaps I should check if it has already been added
            boundaryNodes.Add(node);
            allNodes.Add(node);
        }

        /// <summary>
        /// Checks if the element is internal to this subdomain. If it is, it is added and true is returned. Otherwise false is
        /// returned.
        /// </summary>
        /// <param name="element"></param>
        /// <returns></returns>
        public bool AddElementIfInternal(XContinuumElement2D element)
        {
            //TODO: check if one node is internal, while another is external.
            foreach (var node in element.Nodes)
            {
                if (!(internalNodes.Contains(node) || boundaryNodes.Contains(node))) return false;
            }
            elements.Add(element);
            return true;
        }

        public void AddInternalNode(XNode node)
        {  //TODO: perhaps I should check if it has already been added
            internalNodes.Add(node);
            allNodes.Add(node);
        }

        public int CompareTo(XSubdomain2D_old other)
        {
            return this.ID - other.ID;
        }

        public bool HasEnrichedNodes()
        {
            foreach (var node in allNodes)
            {
                if (node.IsEnriched) return true;
            }
            return false;
        }
    }
}
