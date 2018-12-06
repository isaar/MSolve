using System;
using ISAAR.MSolve.LinearAlgebra.Reduction;

//TODO: perhaps I should return IVectorView instead of IVector. By returning IVectorView I can have classes that only implement
//      IVectorView. On the other hand, I cannot mutate the returned type, so its usefulness is limited.
//TODO: Should IVector Copy() be defined in IVectorView? It doesn't mutate the original vector. However, the return type should
//      still be IVector, since there wouldn't be a point in getting an immutable copy of an immutable class.
namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    /// <summary>
    /// It supports common operations that do not mutate the underlying vector. If you need to store a vector and then pass it
    /// around or allow acceess to it, consider using this interface instead of <see cref="Vector"/> for extra safety.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IVectorView: IIndexable1D, IReducible
    {
        /// <summary>
        /// Performs the following operation for all i:
        /// result[i] = <paramref name="otherCoefficient"/> * <paramref name="otherVector"/>[i] + this[i]. 
        /// Optimized version of <see cref="IVectorView.DoEntrywise(IVectorView, Func{double, double, double})"/> and 
        /// <see cref="IVectorView.LinearCombination(double, IVectorView, double)"/>. Named after BLAS axpy (y = a*x plus y).
        /// The resulting vector is written in a new object and then returned.
        /// </summary>
        /// <param name="otherVector">A vector with the same <see cref="IIndexable1D.Length"/> as this.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if <paramref name="otherVector"/> has different 
        ///     <see cref="IIndexable1D.Length"/> than this.</exception>
        IVector Axpy(IVectorView otherVector, double otherCoefficient);

        /// <summary>
        /// Copies this <see cref="IVector"/> object. A new vector of the same type as this object is initialized and returned.
        /// </summary>
        /// <param name="copyIndexingData">If true, all data of this object will be copied. If false, only the array(s) 
        ///     containing the values of the stored vector entries will be copied. The new vector will reference the same 
        ///     indexing arrays as this one.</param>
        IVector Copy(bool copyIndexingData = false);

        /// <summary>
        /// Initializes a new instance of the same type as this vector, with the exact same storage format and zero entries.
        /// </summary>
        IVector CreateZeroVectorWithSameFormat();

        /// <summary>
        /// Returns an array with the entries of the vector. This is a deep copy operation. For vectors that do not explicitly 
        /// store zeros, the copy may be much larger than the original vector. 
        /// </summary>
        double[] CopyToArray();

        /// <summary>
        /// Performs a binary operation on each pair of entries: 
        /// result[i] = <paramref name="binaryOperation"/>(this[i], <paramref name="vector"/>[i]). 
        /// The resulting vector is written in a new object and then returned.
        /// </summary>
        /// <param name="vector">A vector with the same <see cref="IIndexable1D.Length"/> as this.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if <paramref name="vector"/> has different 
        ///     <see cref="IIndexable1D.Length"/> than this.</exception>
        IVector DoEntrywise(IVectorView vector, Func<double, double, double> binaryOperation);

        /// <summary>
        /// Performs a unary operation on each entry: result[i] = <paramref name="unaryOperation"/>(this[i]).
        /// The resulting vector is written in a new object and then returned.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        IVector DoToAllEntries(Func<double, double> unaryOperation);
        
        /// <summary>
        /// Calculates the dot (or inner/scalar) product of this vector with <paramref name="vector"/>:
        /// result = sum over all i of this[i] * <paramref name="vector"/>[i]).
        /// </summary>
        /// <param name="vector">A vector with the same <see cref="IIndexable1D.Length"/> as this.</param>
        double DotProduct(IVectorView vector);

        /// <summary>
        /// Performs the following operation for all i:
        /// result[i] = <paramref name="thisCoefficient"/> * this[i] + <paramref name="otherCoefficient"/> * 
        /// <paramref name="otherVector"/>[i].
        /// Optimized version of <see cref="DoEntrywiseIntoThis(IVectorView, Func{double, double, double})"/>.
        /// The resulting vector is written in a new object and then returned.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this vector.</param>
        /// <param name="otherVector">A vector with the same <see cref="IIndexable1D.Length"/> as this.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if <paramref name="otherVector"/> has different 
        ///     <see cref="IIndexable1D.Length"/> than this.</exception>
        IVector LinearCombination(double thisCoefficient, IVectorView otherVector, double otherCoefficient);

        /// <summary>
        /// Performs the following operation for all i: result[i] = <paramref name="scalar"/> * this[i].
        /// The resulting vector is written in a new object and then returned.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this vector.</param>
        IVector Scale(double scalar);
    }
}
