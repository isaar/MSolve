using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Enrichments.Tip
{
    interface ITipEnrichment
    {
        double ValueAt(double r, double a);
        double RadialDerivativeAt(double r, double a);
        double AngularDerivativeAt(double r, double a);
    }
}
