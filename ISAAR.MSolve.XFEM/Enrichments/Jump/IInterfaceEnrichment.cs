using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Enrichments.Jump
{
    interface IInterfaceEnrichment
    {
        double ValueAt(double signedDistance);
        double DerivativeAt(double signedDistance);
    }
}
