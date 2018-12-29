//TODO: Move the operators here when C# supports extension operators
using System;

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

        /// <summary>
        /// Performs the following operation for <paramref name="length"/> consecutive entries starting from the provided 
        /// indices: this[i] = this[i] + <paramref name="sourceVector"/>[i].
        /// </summary>
        /// <param name="destinationIndex">The index into this <see cref="IVector"/> where to start overwritting. Constraints:
        ///     <paramref name="destinationIndex"/> + <paramref name="length"/> &lt;= this.<see cref="IIndexable1D.Length"/>.
        ///     </param>
        /// <param name="sourceVector">The other vector operand.</param>
        /// <param name="sourceIndex">The index into <paramref name="sourceVector"/> where to start operating. Constraints: 
        ///     <paramref name="sourceIndex"/> + <paramref name="length"/> &lt;= 
        ///     <paramref name="sourceVector"/>.<see cref="IIndexable1D.Length"/>.</param>
        /// <param name="length">The number of entries to copy.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if <paramref name="length"/> and 
        ///     <paramref name="destinationIndex"/> or <paramref name="sourceIndex"/> violate the described constraints.
        ///     </exception>
        /// <exception cref="Exceptions.PatternModifiedException">Thrown if an entry this[i] needs to be overwritten, but that 
        ///     is not permitted by the vector storage format.</exception>
        public static void AddSubvectorIntoThis(this IVector destinationVector, int destinationIndex, IVectorView sourceVector,
            int sourceIndex, int length) 
            => destinationVector.AxpySubvectorIntoThis(destinationIndex, sourceVector, 1.0, sourceIndex, length);

        // TODO: implement this in each concrete vector.
        public static double Norm2(this IVectorView thisVector) => Math.Sqrt(thisVector.DotProduct(thisVector));

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

        /// <summary>
        /// Performs the following operation for <paramref name="length"/> consecutive entries starting from the provided 
        /// indices: this[i] = this[i] - <paramref name="sourceVector"/>[i].
        /// </summary>
        /// <param name="destinationIndex">The index into this <see cref="IVector"/> where to start overwritting. Constraints:
        ///     <paramref name="destinationIndex"/> + <paramref name="length"/> &lt;= this.<see cref="IIndexable1D.Length"/>.
        ///     </param>
        /// <param name="sourceVector">The other vector operand.</param>
        /// <param name="sourceIndex">The index into <paramref name="sourceVector"/> where to start operating. Constraints: 
        ///     <paramref name="sourceIndex"/> + <paramref name="length"/> &lt;= 
        ///     <paramref name="sourceVector"/>.<see cref="IIndexable1D.Length"/>.</param>
        /// <param name="length">The number of entries to copy.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if <paramref name="length"/> and 
        ///     <paramref name="destinationIndex"/> or <paramref name="sourceIndex"/> violate the described constraints.
        ///     </exception>
        /// <exception cref="Exceptions.PatternModifiedException">Thrown if an entry this[i] needs to be overwritten, but that 
        ///     is not permitted by the vector storage format.</exception>
        public static void SubtractSubvectorIntoThis(this IVector destinationVector, int destinationIndex,
            IVectorView sourceVector, int sourceIndex, int length)
            => destinationVector.AxpySubvectorIntoThis(destinationIndex, sourceVector, -1.0, sourceIndex, length);

        //TODO: remove this
        public static Numerical.LinearAlgebra.Vector ToLegacyVector(this IVectorView vector)
            => Vector.CreateFromVector(vector).ToLegacyVector();

        //TODO: remove this. Its only purpose is to avoid calling IVectorView.CopyToArray() and needless copying, during the
        //      transition phase from the legacy linear algebra design.
        public static double[] ToRawArray(this Vector vector) => vector.RawData;
    }
}
