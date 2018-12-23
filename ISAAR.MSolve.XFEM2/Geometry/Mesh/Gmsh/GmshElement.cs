using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Geometry.Mesh.Gmsh
{
    class GmshElement
    {
        public IsoparametricElementType2D ElementType { get; }
        public IReadOnlyList<XNode2D> Nodes { get; }

        public GmshElement(IsoparametricElementType2D type, IReadOnlyList<XNode2D> nodes)
        {
            this.ElementType = type;
            this.Nodes = nodes;
        }
    }
}
