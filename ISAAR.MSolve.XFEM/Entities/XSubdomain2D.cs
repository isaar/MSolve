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
        private readonly SortedSet<XNode2D> nodes;

        public XSubdomain2D()
        {
            this.nodes = new SortedSet<XNode2D>();
            this.elements = new HashSet<XContinuumElement2D>();
        }

        public ISet<XContinuumElement2D> Elements { get { return elements; } }
        public ISet<XNode2D> Nodes { get; }

        public XSubdomainDofOrderer DofOrder { get; set; }

        public void AddNode(XNode2D node) //TODO: perhaps I should only add elements
        {
            bool alreadyExists = nodes.Add(node);
            if (alreadyExists) throw new ArgumentException("There is already a node with id = " + node.ID);
        }

        public void AddElement(XContinuumElement2D element)
        {
            bool alreadyExists = elements.Add(element);
            if (alreadyExists) throw new ArgumentException("This element is already inserted");
        }

        public void EnumerateDofs()
        {
            DofOrder = XSubdomainDofOrderer.CreateNodeMajor(nodes);
        }
    }
}
