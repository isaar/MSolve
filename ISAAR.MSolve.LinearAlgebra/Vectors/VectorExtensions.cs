//TODO: Move the operators here when C# supports extension operators
namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    /// <summary>
    /// Defines common vector operation shortcuts that can be used as extensions for <see cref="Vector"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class VectorExtensions
    {
        /// <summary>
        /// Performs the operation: result[i] = <paramref name="vector1"/>[i] + <paramref name="vector2"/>[i], 
        /// for 0 &lt;= i &lt; <paramref name="vector1"/>.<see cref="Length"/> = <paramref name="vector2"/>.<see cref="Length"/>.
        /// The resulting vector overwrites the entries of this <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="vector1">The first <see cref="Vector"/> operand. It must have the same <see cref="Length"/> as 
        ///     <paramref name="vector2"/>.</param>
        /// <param name="vector2">The second <see cref="Vector"/> operand. It must have the same <see cref="Length"/> as 
        ///     <paramref name="vector1"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="vector1"/> and <paramref name="vector2"/>
        ///     have different <see cref="Length"/>.</exception>
        public static void AddIntoThis(this Vector vector1, Vector vector2) => vector1.AxpyIntoThis(vector2, 1.0);

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="vector1"/>[i] - <paramref name="vector2"/>[i], 
        /// for 0 &lt;= i &lt; <paramref name="vector1"/>.<see cref="Length"/> = <paramref name="vector2"/>.<see cref="Length"/>.
        /// The resulting vector overwrites the entries of this <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="vector1">The first <see cref="Vector"/> operand. It must have the same <see cref="Length"/> as 
        ///     <paramref name="vector2"/>.</param>
        /// <param name="vector2">The second <see cref="Vector"/> operand. It must have the same <see cref="Length"/> as 
        ///     <paramref name="vector1"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="vector1"/> and <paramref name="vector2"/>
        ///     have different <see cref="Length"/>.</exception>
        public static void SubtractIntoThis(this Vector vector1, Vector vector2) => vector1.AxpyIntoThis(vector2, -1.0);
    }
}
