using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Interpolation;

namespace ISAAR.MSolve.XFEM.Elements
{
    interface IFiniteElement2DView<TNode, TInterpolation>
    {
        IReadOnlyList<TNode> Nodes { get; }
        TInterpolation Interpolation { get; }
    }
}
