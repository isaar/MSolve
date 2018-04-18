using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Reduction
{
    public static class Reductions
    {
        public static double Average(this IVectorView vector)
        {
            return vector.Reduce(0.0, (x, sum) => x + sum, (nz, sum) => sum, sum => sum / vector.Length);
        }

        public static double Max(this IReducible reducible)
        {
            return reducible.Reduce(double.MinValue, (x, max) => x > max ? x : max,
                (nz, max) => 0.0 > max ? 0.0 : max, max => max);
        }

        public static double Min(this IReducible reducible)
        {
            return reducible.Reduce(double.MaxValue, (x, min) => x < min ? x : min,
                (nz, min) => 0.0 < min ? 0.0 : min, min => min);
        }

        // Obsolete? The MKL version should be prefered for performance.
        public static double Norm2(this IReducible reducible)
        {
            return reducible.Reduce(0.0, (x, sum) => x + sum, (nz, sum) => sum, sum => Math.Sqrt(sum));
        }

        public static double Product(this IReducible reducible)
        {
            return reducible.Reduce(1.0, (x, prod) => x * prod, (nz, prod) => nz > 0 ? 0 : prod, prod => prod);
        }

        public static double Sum(this IReducible reducible)
        {
            return reducible.Reduce(0.0, (x, sum) => x + sum, (nz, sum) => sum, sum => sum);
        }
    }
}
