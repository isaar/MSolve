using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Integration.ShapeFunctions
{
    interface IShapeFunction2D
    {
        double ValueAt(double x, double eta);
        double XiDerivativeAt(double x, double y);
        double EtaDerivativeAt(double x, double y);
    }
}
