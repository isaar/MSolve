using System;

//TODO: perhaps some of these methods should only belong to a few concrete classes.
//TODO: perhaps the NonContiguously where one of the two vectors is contiguous, must be bidirectional by setting a flag.
namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    /// <summary>
    /// Operations specified by this interface modify the vector. Therefore it is possible that they may throw exceptions if they 
    /// are used on sparse vector formats and the zero entries are overwritten.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IVector: IVectorView
    {
        /// <summary>
        /// Adds selected entries from <paramref name="otherVector"/> to this vector:
        /// this[<paramref name="thisIndices"/>[i]] += <paramref name="otherVector"/>[<paramref name="otherIndices"/>[i]],
        /// for 0 &lt;= i &lt; <paramref name="thisIndices"/>.Length = <paramref name="otherIndices"/>.Length.
        /// </summary>
        /// <param name="thisIndices">
        /// The indices of this vector, where entries will be added to. Constraints: 
        /// 1) <paramref name="thisIndices"/>.Length == <paramref name="otherIndices"/>.Length,
        /// 2) 0 &lt;= <paramref name="thisIndices"/>[i] &lt; this.<see cref="IIndexable1D.Length"/>, for all valid i
        /// </param>
        /// <param name="otherVector">The vector from which entries will be added.</param>
        /// <param name="otherIndices">
        /// The indices of <paramref name="otherVector"/>, from where entries will be added. Constraints: 
        /// 1) <paramref name="otherIndices"/>.Length == <paramref name="thisIndices"/>.Length,
        /// 2) 0 &lt;= <paramref name="otherIndices"/>[i] &lt; otherVector.<see cref="IIndexable1D.Length"/>, for all valid i
        /// </param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="thisIndices"/> or <paramref name="otherIndices"/> violate the described constraints.
        /// </exception>
        /// <exception cref="IndexOutOfRangeException">
        /// Thrown if <paramref name="thisIndices"/> or <paramref name="otherIndices"/> violate the described constraints.
        /// </exception>
        void AddNonContiguouslyFrom(int[] thisIndices, IVectorView otherVector, int[] otherIndices);

        /// <summary>
        /// Adds selected entries from <paramref name="otherVector"/> to this vector:
        /// this[<paramref name="thisIndices"/>[i]] += <paramref name="otherVector"/>[i], for 0 &lt;= i
        /// &lt; <paramref name="otherVector"/>.<see cref="IIndexable1D.Length"/> = <paramref name="thisIndices"/>.Length.
        /// Contrary to <see cref="AddNonContiguouslyFrom(int[], IVectorView, int[])"/>, access to the entries of 
        /// <paramref name="otherVector"/> is contiguous.
        /// </summary>
        /// <param name="thisIndices">
        /// The indices of this vector, where entries will be added to. Constraints: 
        /// 1) <paramref name="thisIndices"/>.Length == <paramref name="otherVector"/>.<see cref="IIndexable1D.Length"/>,
        /// 2) 0 &lt;= <paramref name="thisIndices"/>[i] &lt; this.<see cref="IIndexable1D.Length"/>, for all valid i
        /// </param>
        /// <param name="otherVector">The vector from which entries will be added.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="thisIndices"/> violates the described constraints.
        /// </exception>
        /// <exception cref="IndexOutOfRangeException">
        /// Thrown if <paramref name="thisIndices"/> violates the described constraints.
        /// </exception>
        void AddNonContiguouslyFrom(int[] thisIndices, IVectorView otherVector);

        /// <summary>
        /// Performs the following operation for all i:
        /// this[i] = <paramref name="otherCoefficient"/> * <paramref name="otherVector"/>[i] + this[i]. 
        /// Optimized version of <see cref="IVector.DoEntrywise(IVectorView, Func{double, double, double})"/> and 
        /// <see cref="IVector.LinearCombination(double, IVectorView, double)"/>. Named after BLAS axpy (y = a*x plus y).
        /// The resulting vector overwrites the entries of this.
        /// </summary>
        /// <param name="otherVector">A vector with the same <see cref="IIndexable1D.Length"/> as this.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="otherVector"/> has different <see cref="IIndexable1D.Length"/> than this.
        /// </exception>
        /// <exception cref="Exceptions.PatternModifiedException">
        /// Thrown if an entry this[i] needs to be overwritten, but that is not permitted by the vector storage format.
        /// </exception> 
        void AxpyIntoThis(IVectorView otherVector, double otherCoefficient);

        /// <summary>
        /// Performs the following operation for <paramref name="length"/> consecutive entries starting from the provided 
        /// indices: this[i] = <paramref name="sourceCoefficient"/> * <paramref name="sourceVector"/>[i] + this[i].
        /// </summary>
        /// <param name="destinationIndex">
        /// The index into this <see cref="IVector"/> where to start overwritting. Constraints:
        /// <paramref name="destinationIndex"/> + <paramref name="length"/> &lt;= this.<see cref="IIndexable1D.Length"/>.
        /// </param>
        /// <param name="sourceVector">The other vector operand.</param>
        /// <param name="sourceCoefficient">A scalar that multiplies each entry of <paramref name="sourceVector"/>.</param>
        /// <param name="sourceIndex">
        /// The index into <paramref name="sourceVector"/> where to start operating. 
        /// Constraints: <paramref name="sourceIndex"/> + <paramref name="length"/> 
        /// &lt;= <paramref name="sourceVector"/>.<see cref="IIndexable1D.Length"/>.
        /// </param>
        /// <param name="length">The number of entries to copy.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="length"/> and <paramref name="destinationIndex"/> or <paramref name="sourceIndex"/> 
        /// violate the described constraints.
        /// </exception>
        /// <exception cref="Exceptions.PatternModifiedException">
        /// Thrown if an entry this[i] needs to be overwritten, but that is not permitted by the vector storage format.
        /// </exception>
        void AxpySubvectorIntoThis(int destinationIndex, IVectorView sourceVector, double sourceCoefficient, int sourceIndex,
            int length);

        /// <summary>
        /// Sets all entries to 0. For sparse or block vectors: the indexing arrays will not be mutated. Therefore the sparsity  
        /// pattern will be preserved. The non-zero entries will be set to 0, but they will still be stored explicitly. 
        /// </summary>
        void Clear();

        /// <summary>
        /// Copies all entries from <paramref name="sourceVector"/> to this <see cref="IVector"/>.
        /// </summary>
        /// <param name="sourceVector">The vector containing the entries to be copied.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="sourceVector"/> has different <see cref="IIndexable1D.Length"/> than this.
        /// </exception>
        /// <exception cref="Exceptions.PatternModifiedException">
        /// Thrown if an entry this[i] needs to be overwritten, but that is not permitted by the vector storage format.
        /// </exception>
        void CopyFrom(IVectorView sourceVector);

        /// <summary>
        /// Copies selected entries from <paramref name="otherVector"/> to this vector:
        /// this[<paramref name="thisIndices"/>[i]] = <paramref name="otherVector"/>[<paramref name="otherIndices"/>[i]],
        /// for 0 &lt;= i &lt; <paramref name="thisIndices"/>.Length = <paramref name="otherIndices"/>.Length.
        /// </summary>
        /// <param name="thisIndices">
        /// The indices of this vector, where entries will be copied to. Constraints: 
        /// 1) <paramref name="thisIndices"/>.Length == <paramref name="otherIndices"/>.Length,
        /// 2) 0 &lt;= <paramref name="thisIndices"/>[i] &lt; this.<see cref="IIndexable1D.Length"/>, for all valid i
        /// </param>
        /// <param name="otherVector">The vector from which entries will be copied.</param>
        /// <param name="otherIndices">
        /// The indices of <paramref name="otherVector"/>, from where entries will be copied. Constraints: 
        /// 1) <paramref name="otherIndices"/>.Length == <paramref name="thisIndices"/>.Length,
        /// 2) 0 &lt;= <paramref name="otherIndices"/>[i] &lt; otherVector.<see cref="IIndexable1D.Length"/>, for all valid i
        /// </param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="thisIndices"/> or <paramref name="otherIndices"/> violate the described constraints.
        /// </exception>
        /// <exception cref="IndexOutOfRangeException">
        /// Thrown if <paramref name="thisIndices"/> or <paramref name="otherIndices"/> violate the described constraints.
        /// </exception>
        void CopyNonContiguouslyFrom(int[] thisIndices, IVectorView otherVector, int[] otherIndices);

        /// <summary>
        /// Copies selected entries from <paramref name="otherVector"/> to this vector:
        /// this[i] = <paramref name="otherVector"/>[<paramref name="otherIndices"/>[i]],
        /// for 0 &lt;= i &lt; this.<see cref="IIndexable1D.Length"/> = <paramref name="otherIndices"/>.Length.
        /// Contrary to <see cref="CopyNonContiguouslyFrom(int[], IVectorView, int[])"/>, access to the entries of this vector
        /// is contiguous.
        /// </summary>
        /// <param name="otherVector">The vector from which entries will be copied.</param>
        /// <param name="otherIndices">
        /// The indices of <paramref name="otherVector"/>, from where entries will be copied. Constraints: 
        /// 1) <paramref name="otherIndices"/>.Length == this.<see cref="IIndexable1D.Length"/>,
        /// 2) 0 &lt;= <paramref name="otherIndices"/>[i] &lt; otherVector.<see cref="IIndexable1D.Length"/>, for all valid i
        /// </param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="otherIndices"/> violates the described constraints.
        /// </exception>
        /// <exception cref="IndexOutOfRangeException">
        /// Thrown if <paramref name="otherIndices"/> violates the described constraints.
        /// </exception>
        void CopyNonContiguouslyFrom(IVectorView otherVector, int[] otherIndices);

        /// <summary>
        /// Copies <paramref name="length"/> consecutive entries from <paramref name="sourceVector"/> to this 
        /// <see cref="IVector"/> starting from the provided indices.
        /// </summary>
        /// <param name="destinationIndex">
        /// The index into this <see cref="IVector"/> where to start copying to. Constraints:
        /// <paramref name="destinationIndex"/> + <paramref name="length"/> &lt;= this.<see cref="IIndexable1D.Length"/>.
        /// </param>
        /// <param name="sourceVector">The vector containing the entries to be copied.</param>
        /// <param name="sourceIndex">
        /// The index into <paramref name="sourceVector"/> where to start copying from. 
        /// Constraints: <paramref name="sourceIndex"/> + <paramref name="length"/> 
        /// &lt;= <paramref name="sourceVector"/>.<see cref="IIndexable1D.Length"/>.
        /// </param>
        /// <param name="length">The number of entries to copy.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="length"/> and <paramref name="destinationIndex"/> or <paramref name="sourceIndex"/> 
        /// violate the described constraints.
        /// </exception>
        /// <exception cref="Exceptions.PatternModifiedException">
        /// Thrown if an entry this[i] needs to be overwritten, but that is not permitted by the vector storage format.
        /// </exception>
        void CopySubvectorFrom(int destinationIndex, IVectorView sourceVector, int sourceIndex, int length);

        /// <summary>
        /// Performs a binary operation on each pair of entries: 
        /// this[i] = <paramref name="binaryOperation"/>(this[i], <paramref name="vector"/>[i]). 
        /// The resulting vector overwrites the entries of this.
        /// </summary>
        /// <param name="vector">A vector with the same <see cref="IIndexable1D.Length"/> as this.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="vector"/> has different <see cref="IIndexable1D.Length"/> than this.
        /// </exception>
        /// <exception cref="Exceptions.PatternModifiedException">
        /// Thrown if an entry this[i] needs to be overwritten, but that is not permitted by the vector storage format.
        /// </exception> 
        void DoEntrywiseIntoThis(IVectorView otherVector, Func<double, double, double> binaryOperation);

        /// <summary>
        /// Performs a unary operation on each entry: this[i] = <paramref name="unaryOperation"/>(this[i]).
        /// The resulting vector overwrites the entries of this.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        /// <exception cref="Exceptions.PatternModifiedException">
        /// Thrown if an entry this[i] needs to be overwritten, but that is not permitted by the vector storage format.
        /// </exception> 
        void DoToAllEntriesIntoThis(Func<double, double> unaryOperation);

        /// <summary>
        /// Performs the following operation for all i:
        /// this[i] = <paramref name="thisCoefficient"/> * this[i] + <paramref name="otherCoefficient"/> * 
        /// <paramref name="otherMatrix"/>[i].
        /// Optimized version of <see cref="DoEntrywiseIntoThis(IVectorView, Func{double, double, double})"/>.
        /// The resulting vector overwrites the entries of this.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this vector.</param>
        /// <param name="otherVector">A vector with the same <see cref="IIndexable1D.Length"/> as this.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="otherVector"/> has different <see cref="IIndexable1D.Length"/> than this.
        /// </exception>
        /// <exception cref="Exceptions.PatternModifiedException">
        /// Thrown if an entry this[i] needs to be overwritten, but that is not permitted by the vector storage format.
        /// </exception> 
        void LinearCombinationIntoThis(double thisCoefficient, IVectorView otherVector, double otherCoefficient);

        /// <summary>
        /// Performs the following operation for all i: this[i] = <paramref name="scalar"/> * this[i].
        /// The resulting vector overwrites the entries of this.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this vector.</param>
        void ScaleIntoThis(double scalar);

        /// <summary>
        /// Setter that will work as expected for general dense vectors. For sparse vectors it will throw a 
        /// <see cref="Exceptions.SparsityPatternModifiedException"/>, if a structural zero entry is written to.
        /// </summary>
        /// <param name="index">
        /// The index of the entry to set. Constraints: 0 &lt;= <paramref name="index"/> &lt; <see cref="IIndexable1D.Length"/>.
        /// </param>
        /// <param name="value">The new value of this[<paramref name="index"/>].</param>
        /// <exception cref="IndexOutOfRangeException">
        /// Thrown if <paramref name="index"/> violates the described constraints.
        /// </exception>
        /// <exception cref="Exceptions.SparsityPatternModifiedException"> 
        /// Thrown if a structural zero entry of a sparse vector format is written to.
        /// </exception>
        void Set(int index, double value);
    }
}
