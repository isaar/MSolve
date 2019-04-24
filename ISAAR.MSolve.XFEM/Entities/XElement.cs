using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.XFEM.Elements;

namespace ISAAR.MSolve.XFEM.Entities
{
    public class XElement : IElement
    {
        private readonly int id;

        public XElement(int id, IXFiniteElement elementType)
        {
            // TODO: Add checks here
            this.id = id;
            this.ElementType = elementType;
        }

        public int ID
        {
            get => id;
            set => throw new InvalidOperationException("The id of an element cannot be changed");
        }

        IElementType IElement.ElementType => ElementType;
        public IXFiniteElement ElementType { get; }

        IReadOnlyList<INode> IElement.Nodes => ElementType.Nodes;
        public IReadOnlyList<XNode> Nodes => ElementType.Nodes;

        ISubdomain IElement.Subdomain => this.Subdomain;
        public XSubdomain Subdomain { get; set; }
    }
}
