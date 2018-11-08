using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Elements
{
    interface IFiniteElement2DView<TNode, TInterpolation>
    {
        IReadOnlyList<TNode> Nodes { get; }
        TInterpolation Interpolation { get; }
    }
}
