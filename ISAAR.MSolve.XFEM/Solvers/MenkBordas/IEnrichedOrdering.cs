using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    interface IEnrichedOrdering
    {
        void ReorderEnrichedDofs(XSubdomain2D subdomainDofOrderer);
    }
}
