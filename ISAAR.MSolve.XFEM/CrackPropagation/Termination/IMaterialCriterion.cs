using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.CrackPropagation.Termination
{
    interface IMaterialCriterion
    {
        bool Terminate(double sifMode1, double sifMode2);
    }
}
