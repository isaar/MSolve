using System;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Perhaps I should delete Norm2() from here. It is implemented by concrete classes much more efficiently.
namespace ISAAR.MSolve.LinearAlgebra.Reduction
{
    /// <summary>
    /// Defines a number of common reduction operations that can be used as extensions for <see cref="IReducible"/>. For more 
    /// see https://en.wikipedia.org/wiki/Fold_(higher-order_function).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class Reductions
    {
        /// <summary>
        /// Calculates the average over all entries of a vector.
        /// </summary>
        /// <param name="vector">The vector upon which this reduction will be applied.</param>
        public static double Average(this IVectorView vector)
        {
            return vector.Reduce(0.0, (x, sum) => x + sum, (nz, sum) => sum, sum => sum / vector.Length);
        }

        /// <summary>
        /// Calculates the maximum value of all entries of an <see cref="IReducible"/>.
        /// </summary>
        /// <param name="reducible">A matrix, vector or similar collection.</param>
        public static double Max(this IReducible reducible)
        {
            return reducible.Reduce(double.MinValue, (x, max) => x > max ? x : max,
                (nz, max) => 0.0 > max ? 0.0 : max, max => max);
        }

        /// <summary>
        /// Calculates the maximum absolute value of all entries of an <see cref="IReducible"/>.
        /// </summary>
        /// <param name="reducible">A matrix, vector or similar collection.</param>
        public static double MaxAbsolute(this IReducible reducible)
        {
            return reducible.Reduce(double.MaxValue, (x, max) =>
            {
                double abs = Math.Abs(x);
                return abs > max ? abs : max;
            },
            (nz, max) => 0.0 > max ? 0.0 : max, max => max);
        }

        /// <summary>
        /// Calculates the minimum value of all entries of an <see cref="IReducible"/>.
        /// </summary>
        /// <param name="reducible">A matrix, vector or similar collection.</param>
        public static double Min(this IReducible reducible)
        {
            return reducible.Reduce(double.MaxValue, (x, min) => x < min ? x : min,
                (nz, min) => 0.0 < min ? 0.0 : min, min => min);
        }

        /// <summary>
        /// Calculates the minimum absolute value of all entries of an <see cref="IReducible"/>.
        /// </summary>
        /// <param name="reducible">A matrix, vector or similar collection.</param>
        public static double MinAbsolute(this IReducible reducible)
        {
            return reducible.Reduce(double.MaxValue, (x, min) =>
            {
                double abs = Math.Abs(x);
                return abs < min ? abs : min;
            },
            (nz, min) => 0.0 < min ? 0.0 : min, min => min);
        }

        /// <summary>
        /// Calculates the Euclidian norm or 2-norm over all entries of an <see cref="IReducible"/>. For more see
        /// https://en.wikipedia.org/wiki/Norm_(mathematics)#Euclidean_norm. WARNING: most classes that implement 
        /// <see cref="IReducible"/> will provide far more efficient methods for calculating norms. Use them instead.
        /// </summary>
        /// <param name="reducible">A matrix, vector or similar collection.</param>
        public static double Norm2(this IReducible reducible)
        {
            return reducible.Reduce(0.0, (x, sum) => x + sum, (nz, sum) => sum, sum => Math.Sqrt(sum));
        }

        /// <summary>
        /// Calculates the product over all entries of an <see cref="IReducible"/>.
        /// </summary>
        /// <param name="reducible">A matrix, vector or similar collection.</param>
        public static double Product(this IReducible reducible)
        {
            return reducible.Reduce(1.0, (x, prod) => x * prod, (nz, prod) => nz > 0 ? 0 : prod, prod => prod);
        }

        /// <summary>
        /// Calculates the sum over all entries of an <see cref="IReducible"/>.
        /// </summary>
        /// <param name="reducible">A matrix, vector or similar collection.</param>
        public static double Sum(this IReducible reducible)
        {
            return reducible.Reduce(0.0, (x, sum) => x + sum, (nz, sum) => sum, sum => sum);
        }
    }
}
