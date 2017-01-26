using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;

namespace ISAAR.MSolve.XFEM.Entities
{
    class Element2D
    {
        public IReadOnlyList<Node2D> Nodes { get; }
        public IFiniteElement2D ElementType { get; private set; }

        public Element2D(IEnumerable<Node2D> nodes, IFiniteElement2D elementType)
        {
            // TODO: Add checks here
            this.Nodes = new List<Node2D>(nodes);
            this.ElementType = elementType;
        }
    }
}
