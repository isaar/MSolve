using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Output.VTK
{
    interface IXfemOutput
    {
        void WriteOutputData(IDofOrderer dofOrderer, Vector freeDisplacements, Vector constrainedDisplacements, int step);
    }
}
