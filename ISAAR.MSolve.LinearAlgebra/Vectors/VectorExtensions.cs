//TODO: Move the operators here when C# supports extension operators
using System;
using ISAAR.MSolve.LinearAlgebra.Exceptions;

//TODO: Use the generic interfaces IEntrywiseOperable1D, etc (and create some for axpy, linear combo), instead of implementing
//      the extensions for each vector/matrix type and the IVectorView, etc. interfaces. However the extension method should be
//      concise as possible. Having to declare the generic types (see the generic MultiplyEntrywise) is prohibitive, especially
//      if IntelliSense does not suggest the generic extension method. How dids LINQ solve this issue?
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

        /// <summary>
        /// Performs the operation: a x b = { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]},
        /// where a = <paramref name="thisVector"/>, b = <paramref name="otherVector"/> and both have exactly 3 entries..
        /// The result is a vector. Also note that: other.Cross(this) = - this.Cross(other).
        /// </summary>
        /// <param name="thisVector">A vector with three entries.</param>
        /// <param name="otherVector">A vector with three entries.</param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="thisVector"/> or <paramref name="otherVector"/> do not have 
        /// <see cref="IIndexable1D.Length"/> = 3.
        /// </exception>
        public static Vector CrossProduct(this Vector thisVector, Vector otherVector) //TODO: Should this be a member method? It breaks encapsulation as it is
        {
            if (thisVector.Length != 3) throw new NonMatchingDimensionsException(
                $"Vector 1 has length = {thisVector.Length} instead of 3");
            if (otherVector.Length != 3) throw new NonMatchingDimensionsException(
                $"Vector 2 has length = {thisVector.Length} instead of 3");

            double[] a = thisVector.RawData;
            double[] b = otherVector.RawData;

            return Vector.CreateFromArray(new double[]
            {
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0]
            });
        }

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="thisVector"/>[i] * <paramref name="otherVector"/>[i], for all
        /// valid i. The resulting vector is written in a new object and then returned.
        /// </summary>
        /// <param name="thisVector">A vector.</param>
        /// <param name="otherVector">A vector with the same <see cref="IIndexable1D.Length"/> as this vector.</param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="otherVector"/> has different <see cref="IIndexable1D.Length"/> than this vector.
        /// </exception>
        public static TVectorOut MultiplyEntrywise<TVectorIn, TVectorOut>(
            this IEntrywiseOperableView1D<TVectorIn, TVectorOut> thisVector, TVectorIn otherVector)
            where TVectorIn: IVectorView
            where TVectorOut: IVector
            => thisVector.DoEntrywise(otherVector, (x, y) => x * y); //TODO: nice in theory, but passing a lambda to DoEntrywise is less verbose.

        /// <summary>
        /// Performs the operation: 
        /// <paramref name="thisVector"/>[i] = <paramref name="thisVector"/>[i] * <paramref name="otherVector"/>[i], for all
        /// valid i. The resulting vector overwrites the entries of this vector.
        /// </summary>
        /// <param name="thisVector">A vector.</param>
        /// <param name="otherVector">A vector with the same <see cref="IIndexable1D.Length"/> as this vector.</param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="otherVector"/> has different <see cref="IIndexable1D.Length"/> than this vector.
        /// </exception>
        /// <exception cref="PatternModifiedException">
        /// Thrown if an entry this[i] needs to be overwritten, but that is not permitted by the vector storage format.
        /// </exception>
        public static void MultiplyEntrywiseIntoThis<TVectorIn>(
            this IEntrywiseOperable1D<TVectorIn> thisVector, TVectorIn otherVector)
            where TVectorIn: IVectorView
            => thisVector.DoEntrywiseIntoThis(otherVector, (x, y) => x * y); //TODO: nice in theory, but passing a lambda to DoEntrywise() is less verbose.

        /// <summary>
        /// Performs the operation: result[i] = this[i] ^ 0.5 for all valid i. 
        /// The resulting vector is written in a new object and then returned.
        /// </summary>
        public static Vector Sqrt(this Vector vector) => vector.DoToAllEntries(x => Math.Sqrt(x));

        /// <summary>
        /// Performs the operation: this[i] = this[i] ^ 0.5 for all valid i. 
        /// The resulting vector overwrites the entries of this vector.
        /// </summary>
        public static void SqrtIntoThis(this Vector vector) => vector.DoToAllEntriesIntoThis(x => Math.Sqrt(x));

        /// <summary>
        /// Performs the operation: result[i] = this[i] ^ 2 for all valid i. 
        /// The resulting vector is written in a new object and then returned.
        /// </summary>
        public static Vector Square(this Vector vector) => vector.DoToAllEntries(x => x * x);

        /// <summary>
        /// Performs the operation: this[i] = this[i] ^ 2 for all valid i. 
        /// The resulting vector overwrites the entries of this vector.
        /// </summary>
        public static void SquareIntoThis(this Vector vector) => vector.DoToAllEntriesIntoThis(x => x * x);

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
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="length"/> and 
        ///     <paramref name="destinationIndex"/> or <paramref name="sourceIndex"/> violate the described constraints.
        ///     </exception>
        /// <exception cref="PatternModifiedException">Thrown if an entry this[i] needs to be overwritten, but that 
        ///     is not permitted by the vector storage format.</exception>
        public static void SubtractSubvectorIntoThis(this IVector destinationVector, int destinationIndex,
            IVectorView sourceVector, int sourceIndex, int length)
            => destinationVector.AxpySubvectorIntoThis(destinationIndex, sourceVector, -1.0, sourceIndex, length);
    }
}
