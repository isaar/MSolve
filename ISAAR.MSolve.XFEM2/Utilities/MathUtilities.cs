using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Utilities
{
    static class MathUtilities
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
        /// Solves a * x^2 + b * x + c = 0
        /// </summary>
        /// <param name="quadCoeff">Coefficient of quadratic term</param>
        /// <param name="linCoeff">Coefficient of linear term</param>
        /// <param name="constCoeff">Coeffecient of constant term</param>
        /// <returns>The solutions in an array of size 0, 1 or 2</returns>
        public static double[] SolveQuadraticEquation(double quadCoeff, double linCoeff, double constCoeff)
        {
            double discriminant = linCoeff * linCoeff - 4 * quadCoeff * constCoeff;
            if (IsZero(discriminant)) return new double[] { -linCoeff / (2 * quadCoeff) };
            else if (discriminant < 0) return new double[0];
            else
            {
                double sqrtD = Math.Sqrt(discriminant);
                double x1 = (-linCoeff + sqrtD) / (2 * quadCoeff);
                double x2 = (-linCoeff - sqrtD) / (2 * quadCoeff);
                return new double[] { x1, x2 };
            }
        }

        private static bool IsZero(double value)
        {
            return Math.Abs(value) <= TOLERANCE;
        }
    }
}
