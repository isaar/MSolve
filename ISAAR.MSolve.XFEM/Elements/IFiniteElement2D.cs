using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.XFEM.Elements
{
    interface IFiniteElement2D
    {
        SymmetricMatrix2D<double> BuildStiffnessMatrix();
    }
}
