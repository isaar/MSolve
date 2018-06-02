using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Entities.Decomposition
{
    interface IDomainDecomposer
    {
        XCluster2D CreateSubdomains();
        void UpdateSubdomains(XCluster2D cluster);
    }
}
