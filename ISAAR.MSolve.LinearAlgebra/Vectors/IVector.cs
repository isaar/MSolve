using System;

namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    /// <summary>
    /// Operations specified by this interface modify the vector. Therefore it is possible that they may throw exceptions if they 
    /// are used on sparse vector formats.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IVector: IVectorView
    {
        /// <summary>
        /// Performs the following operation for all i:
        /// this[i] = <paramref name="otherCoefficient"/> * <paramref name="otherVector"/>[i] + this[i]. 
        /// Optimized version of <see cref="IVector.DoEntrywise(IVectorView, Func{double, double, double})"/> and 
        /// <see cref="IVector.LinearCombination(double, IVectorView, double)"/>. Named after BLAS axpy (y = a*x plus y).
        /// The resulting vector overwrites the entries of this.
        /// </summary>
        /// <param name="otherVector">A vector with the same <see cref="IIndexable1D.Length"/> as this.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if <paramref name="otherVector"/> has different 
        ///     <see cref="IIndexable1D.Length"/> than this.</exception>
        /// <exception cref="Exceptions.PatternModifiedException">Thrown if an entry this[i] needs to be overwritten, but that 
        ///     is not permitted by the vector storage format.</exception> 
        void AxpyIntoThis(IVectorView otherVector, double otherCoefficient);

        /// <summary>
        /// Performs a binary operation on each pair of entries: 
        /// this[i] = <paramref name="binaryOperation"/>(this[i], <paramref name="vector"/>[i]). 
        /// The resulting vector overwrites the entries of this.
        /// </summary>
        /// <param name="vector">A vector with the same <see cref="IIndexable1D.Length"/> as this.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if <paramref name="vector"/> has different 
        ///     <see cref="IIndexable1D.Length"/> than this.</exception>
        /// <exception cref="Exceptions.PatternModifiedException">Thrown if an entry this[i] needs to be overwritten, but that 
        ///     is not permitted by the vector storage format.</exception> 
        void DoEntrywiseIntoThis(IVectorView other, Func<double, double, double> binaryOperation);

        /// <summary>
        /// Performs a unary operation on each entry: this[i] = <paramref name="unaryOperation"/>(this[i]).
        /// The resulting vector overwrites the entries of this.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        /// <exception cref="Exceptions.PatternModifiedException">Thrown if an entry this[i] needs to be overwritten, but that 
        ///     is not permitted by the vector storage format.</exception> 
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
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if <paramref name="otherVector"/> has different 
        ///     <see cref="IIndexable1D.Length"/> than this.</exception>
        /// <exception cref="Exceptions.PatternModifiedException">Thrown if an entry this[i] needs to be overwritten, but that 
        ///     is not permitted by the vector storage format.</exception> 
        void LinearCombinationIntoThis(double thisCoefficient, IVectorView otherVector, double otherCoefficient);

        /// <summary>
        /// Performs the following operation for all i: this[i] = <paramref name="scalar"/> * this[i].
        /// The resulting vector overwrites the entries of this.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this vector.</param>
        void ScaleIntoThis(double scalar);

        /// <summary>
        /// Setter that will work as expected for general dense vectors. For sparse vectors it will throw a 
        /// <see cref="Exceptions.SparsityPatternModifiedException"/> if a structural zero entry is written to.
        /// </summary>
        /// <param name="index">The index of the entry to set. Constraints: 
        ///     0 &lt;= <paramref name="index"/> &lt; <see cref="IIndexable1D.Length"/>.</param>
        /// <param name="value">The new value of this[<paramref name="index"/>].</param>
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="index"/> violates the described 
        ///     constraints.</exception>
        /// <exception cref="Exceptions.SparsityPatternModifiedException"> Thrown if a structural zero entry of a sparse vector 
        ///     format is written to.</exception>
        void SetEntryRespectingPattern(int index, double value);
    }
}
