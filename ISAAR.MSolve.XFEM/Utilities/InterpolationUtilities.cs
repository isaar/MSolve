using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Utilities
{
    static class InterpolationUtilities
    {
        public static double InterpolateNodalValuesToPoint(double[] nodalValues, double[] interpolationFunctionsAtPoint)
        {
            Debug.Assert(nodalValues.Length == interpolationFunctionsAtPoint.Length);
            double sum = 0.0;
            for (int i = 0; i < nodalValues.Length; ++i)
            {
                sum += nodalValues[i] * interpolationFunctionsAtPoint[i];
            }
            return sum;
        }
    }
}
