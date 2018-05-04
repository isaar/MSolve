using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Entities
{
    /// <summary>
    /// Only for enriched dofs.
    /// </summary>
    class XSubdomain2D
    {
        private readonly HashSet<XContinuumElement2D> elements;
        private readonly SortedSet<XNode2D> allNodes; //TODO: Use a SortedSuperSet custom set class.
        private readonly SortedSet<XNode2D> boundaryNodes;
        private readonly SortedSet<XNode2D> internalNodes;

        public XSubdomain2D()
        {
            this.allNodes = new SortedSet<XNode2D>();
            this.elements = new HashSet<XContinuumElement2D>();
            this.boundaryNodes = new SortedSet<XNode2D>();
            this.internalNodes = new SortedSet<XNode2D>();
        }

        public ISet<XContinuumElement2D> Elements { get { return elements; } }
        public ISet<XNode2D> AllNodes { get { return allNodes; } }
        public ISet<XNode2D> BoundaryNodes { get { return boundaryNodes; } }
        public ISet<XNode2D> InternalNodes { get { return internalNodes; } }
        public XSubdomainDofOrderer DofOrderer { get; set; }

        public void AddBoundaryNode(XNode2D node) //TODO: perhaps I should only add elements
        {
            bool isNew = boundaryNodes.Add(node);
            isNew |= allNodes.Add(node);
            if (!isNew) throw new ArgumentException("There is already a node with id = " + node.ID);
        }

        public void AddInternalNode(XNode2D node) //TODO: perhaps I should only add elements
        {
            bool isNew = internalNodes.Add(node);
            isNew |= allNodes.Add(node);
            if (!isNew) throw new ArgumentException("There is already a node with id = " + node.ID);
        }

        public void AddElement(XContinuumElement2D element)
        {
            bool isNew = elements.Add(element);
            if (!isNew) throw new ArgumentException("This element is already inserted");
        }
    }
}
