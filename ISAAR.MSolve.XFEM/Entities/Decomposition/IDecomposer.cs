using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Entities.Decomposition
{
    interface IDecomposer
    {
        XCluster2D CreateSubdomains();
    }
}
