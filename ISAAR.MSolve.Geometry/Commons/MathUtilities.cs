using System;
using System.Collections.Generic;

namespace ISAAR.MSolve.Geometry.Commons
{
    public static class MathUtilities
    {
        private static readonly double TOLERANCE = 1.0e-6;

        public static int IndexOfMinAbs(IReadOnlyList<double> values)
        {
            double min = double.MaxValue;
            int pos = -1;
            for (int i = 0; i < values.Count; ++i)
            {
                double absDistance = Math.Abs(values[i]);
                if (absDistance < min)
                {
                    min = absDistance;
                    pos = i;
                }
            }
            return pos;
        }

        /// <summary>
        /// Wraps a counter-clockwise angle in (-pi, pi]
        /// <returns></returns>
        public static double WrapAngle(double angle)
        {
            // TODO: (-pi, pi] perhaps is not the best range to work with. 
            // It is convenient in that atan2 returns values there and negative angles are non convex, but the interval should 
            // be closed on the lower bound, to match other formulas. Overall [0, 2pi) seems better overall.

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

        private static bool IsZero(double value) => Math.Abs(value) <= TOLERANCE;
    }
}
