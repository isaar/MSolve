using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Enrichments.Jump
{
    class Ramp: IInterfaceEnrichment
    {
        public double ValueAt(double signedDistance)
        {
            return Math.Abs(signedDistance);
        }

        public double DerivativeAt(double signedDistance)
        {

            return Math.Sign(signedDistance);
        }
    }
}
