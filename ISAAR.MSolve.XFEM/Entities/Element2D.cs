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
        public int ID { get; }
        public IReadOnlyList<XNode2D> Nodes { get { return ElementType.Nodes; } }
        public XElement2D ElementType { get; }

        public Element2D(int id, XElement2D elementType)
        {
            // TODO: Add checks here
            this.ID = id;
            this.ElementType = elementType;
        }
    }
}
