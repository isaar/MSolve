using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    interface IStandardOrdering
    {
        void ReorderStandardDofs(XClusterDofOrderer stdDofOrderer);
    }
}
