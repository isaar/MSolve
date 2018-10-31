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
        /// Performs the operation: 
        /// result[i] = <paramref name="thisVector"/>[i] + <paramref name="otherVector"/>[i], for 
        /// 0 &lt;= i &lt; <paramref name="thisVector"/>.<see cref="IIndexable1D.Length"/> 
        /// = <paramref name="otherVector"/>.<see cref="IIndexable1D.Length"/>.
        /// The resulting vector is written in a new object and then returned.
        /// </summary>
        /// <param name="thisVector">
        /// A vector with the same <see cref="IIndexable1D.Length"/> as <paramref name="otherVector"/>.
        /// </param>
        /// <param name="otherVector">
        /// A vector with the same <see cref="IIndexable1D.Length"/> as <paramref name="thisVector"/>.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="thisVector"/> and <paramref name="otherVector"/> have different 
        /// <see cref="IIndexable1D.Length"/>.
        /// </exception>
        public static IVector Add(this IVectorView thisVector, IVectorView otherVector)
            => thisVector.Axpy(otherVector, 1.0);

        /// <summary>
        /// Performs the operation: 
        /// <paramref name="thisVector"/>[i] = <paramref name="thisVector"/>[i] + <paramref name="otherVector"/>[i], for 
        /// 0 &lt;= i &lt; <paramref name="thisVector"/>.<see cref="IIndexable1D.Length"/> 
        /// = <paramref name="otherVector"/>.<see cref="IIndexable1D.Length"/>.
        /// The resulting vector overwrites the entries of this <see cref="IVector"/> instance.
        /// </summary>
        /// <param name="thisVector">
        /// A vector with the same <see cref="IIndexable1D.Length"/> as <paramref name="otherVector"/>.
        /// </param>
        /// <param name="otherVector">
        /// A vector with the same <see cref="IIndexable1D.Length"/> as <paramref name="thisVector"/>.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="thisVector"/> and <paramref name="otherVector"/> have different 
        /// <see cref="IIndexable1D.Length"/>.
        /// </exception>
        public static void AddIntoThis(this IVector thisVector, IVectorView otherVector) 
            => thisVector.AxpyIntoThis(otherVector, 1.0);


        /// <summary>
        /// Performs the operation: 
        /// <paramref name="thisVector"/>[i] = <paramref name="thisVector"/>[i] + <paramref name="otherVector"/>[i], for 
        /// 0 &lt;= i &lt; <paramref name="thisVector"/>.<see cref="IIndexable1D.Length"/> 
        /// = <paramref name="otherVector"/>.<see cref="IIndexable1D.Length"/>.
        /// The resulting vector overwrites the entries of this <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="thisVector">
        /// A vector with the same <see cref="IIndexable1D.Length"/> as <paramref name="otherVector"/>.
        /// </param>
        /// <param name="otherVector">
        /// A vector with the same <see cref="IIndexable1D.Length"/> as <paramref name="thisVector"/>.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="thisVector"/> and <paramref name="otherVector"/> have different 
        /// <see cref="IIndexable1D.Length"/>.
        /// </exception>
        public static void AddIntoThis(this Vector thisVector, Vector otherVector) 
            => thisVector.AxpyIntoThis(otherVector, 1.0);

        // TODO: implement this in each concrete vector.
        public static double Norm2(this IVectorView thisVector) => thisVector.DotProduct(thisVector);

        /// <summary>
        /// Performs the operation: 
        /// result[i] = <paramref name="thisVector"/>[i] - <paramref name="otherVector"/>[i], for 
        /// 0 &lt;= i &lt; <paramref name="thisVector"/>.<see cref="IIndexable1D.Length"/> 
        /// = <paramref name="otherVector"/>.<see cref="IIndexable1D.Length"/>.
        /// The resulting vector is written in a new object and then returned.
        /// </summary>
        /// <param name="thisVector">
        /// A vector with the same <see cref="IIndexable1D.Length"/> as <paramref name="otherVector"/>.
        /// </param>
        /// <param name="otherVector">
        /// A vector with the same <see cref="IIndexable1D.Length"/> as <paramref name="thisVector"/>.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="thisVector"/> and <paramref name="otherVector"/> have different 
        /// <see cref="IIndexable1D.Length"/>.
        /// </exception>
        public static IVector Subtract(this IVectorView thisVector, IVectorView otherVector)
            => thisVector.Axpy(otherVector, -1.0);

        /// <summary>
        /// Performs the operation: 
        /// <paramref name="thisVector"/>[i] = <paramref name="thisVector"/>[i] - <paramref name="otherVector"/>[i], 
        /// for 0 &lt;= i &lt; <paramref name="thisVector"/>.<see cref="IIndexable1D.Length"/> 
        /// = <paramref name="otherVector"/>.<see cref="IIndexable1D.Length"/>.
        /// The resulting vector overwrites the entries of this <see cref="IVector"/> instance.
        /// </summary>
        /// <param name="thisVector">
        /// A vector with the same <see cref="IIndexable1D.Length"/> as <paramref name="otherVector"/>.
        /// </param>
        /// <param name="otherVector">
        /// A vector with the same <see cref="IIndexable1D.Length"/> as <paramref name="thisVector"/>.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="thisVector"/> and <paramref name="otherVector"/> have different 
        /// <see cref="IIndexable1D.Length"/>.
        /// </exception>
        public static void SubtractIntoThis(this IVector thisVector, IVectorView otherVector)
            => thisVector.AxpyIntoThis(otherVector, -1.0);

        /// <summary>
        /// Performs the operation: 
        /// <paramref name="thisVector"/>[i] = <paramref name="thisVector"/>[i] - <paramref name="otherVector"/>[i], 
        /// for 0 &lt;= i &lt; <paramref name="thisVector"/>.<see cref="IIndexable1D.Length"/> 
        /// = <paramref name="otherVector"/>.<see cref="IIndexable1D.Length"/>.
        /// The resulting vector overwrites the entries of this <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="thisVector">
        /// A vector with the same <see cref="IIndexable1D.Length"/> as <paramref name="otherVector"/>.
        /// </param>
        /// <param name="otherVector">
        /// A vector with the same <see cref="IIndexable1D.Length"/> as <paramref name="thisVector"/>.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="thisVector"/> and <paramref name="otherVector"/> have different 
        /// <see cref="IIndexable1D.Length"/>.
        /// </exception>
        public static void SubtractIntoThis(this Vector thisVector, Vector otherVector) 
            => thisVector.AxpyIntoThis(otherVector, -1.0);

        //TODO: remove this
        public static Numerical.LinearAlgebra.Vector ToLegacyVector(this IVectorView vector)
            => Vector.CreateFromVector(vector).ToLegacyVector();
    }
}
