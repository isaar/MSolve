using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Utilities
{
    // TODO: (-pi, pi] perhaps is not the best range to work with. 
    // It is convenient in that atan2 returns values there and negative angles are non convex, but the interval should 
    // be closed on the lower bound, to match other formulas. Overall [0, 2pi) seems better overall.
    static class AngleUtilities
    {
        /// <summary>
        /// Wraps a counter-clockwise angle in (-pi, pi]
        /// <returns></returns>
        public static double Wrap(double angle)
        {
            double twoPI = 2.0 * Math.PI;
            // Wrap to [0, 2pi)
            double quotient = Math.Floor(angle / twoPI);
            double modulus = angle - twoPI * quotient;
            // Wrap to (-pi, pi]
            double excess = modulus - Math.PI;
            if (excess > 0) // (pi, 2pi) -> (-pi, 0). The [0, pi] is not affected.
            {
                modulus = -Math.PI + excess;
            }
            return modulus;
           
        }
    }
}
