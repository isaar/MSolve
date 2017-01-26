using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Integration
{
    interface IShapeFunction2D
    {
        double ValueAt(double xi, double eta);
        double XiDerivativeAt(double xi, double eta);
        double EtaDerivativeAt(double xi, double eta);
    }
}
