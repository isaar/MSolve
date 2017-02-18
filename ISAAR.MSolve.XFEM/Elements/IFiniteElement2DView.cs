using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Elements
{
    interface IFiniteElement2DView
    {
        IReadOnlyList<Node2D> Nodes { get; }

        SymmetricMatrix2D<double> BuildStiffnessMatrix();
    }
}
