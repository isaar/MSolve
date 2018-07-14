using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//TODO: Move the operators here when C# supports extension operators
namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    public static class VectorExtensions
    {
        /// <summary>
        /// result = vector1 + vector2
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        /// <returns></returns>
        public static Vector Add(this Vector vector1, Vector vector2)
        {
            return vector1.Axpy(1.0, vector2);
        }

        /// <summary>
        /// vector1 = vector1 + vector2
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        public static void AddIntoThis(this Vector vector1, Vector vector2)
        {
            vector1.AxpyIntoThis(vector2, 1.0);
        }

        /// <summary>
        /// Computes the Hadamard product of two vectors: result[i] = this[i] * other[i]. This is not the dot product.
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        /// <returns></returns>
        public static Vector MultiplyPointwise(this Vector vector1, Vector vector2)
        {
            return vector1.DoEntrywise(vector2, (x1, x2) => x1 * x2);
        }

        /// <summary>
        /// Computes the Hadamard product of two vectors: this[i] = this[i] * other[i]. This is not the dot product.
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        public static void MultiplyPointwiseInPlace(this Vector vector1, Vector vector2)
        {
            vector1.DoEntrywiseIntoThis(vector2, (x1, x2) => x1 * x2);
        }

        /// <summary>
        /// result = vector1 - vector1
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        /// <returns></returns>
        public static Vector Subtract(this Vector vector1, Vector vector2)
        {
            return vector1.Axpy(-1.0, vector2);
        }

        /// <summary>
        /// vector1 = vector1 - vector2
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        public static void SubtractIntoThis(this Vector vector1, Vector vector2)
        {
            vector1.AxpyIntoThis(vector2, - 1.0);
        }
    }
}
