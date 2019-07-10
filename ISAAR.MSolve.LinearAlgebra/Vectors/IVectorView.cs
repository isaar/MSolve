using System;
using ISAAR.MSolve.LinearAlgebra.Reduction;

//TODO: perhaps I should return IVectorView instead of IVector. By returning IVectorView I can have classes that only implement
//      IVectorView. On the other hand, I cannot mutate the returned type, so its usefulness is limited.
//TODO: Should IVector Copy() be defined in IVectorView? It doesn't mutate the original vector. However, the return type should
//      still be IVector, since there wouldn't be a point in getting an immutable copy of an immutable class. Same for IMatrixView
namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    /// <summary>
    /// It supports common operations that do not mutate the underlying vector. If you need to store a vector and then pass it
    /// around or allow acceess to it, consider using this interface instead of <see cref="Vector"/> for extra safety.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IVectorView: IIndexable1D, IReducible, IEntrywiseOperableView1D<IVectorView, IVector>
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
        /// Copies this <see cref="IVectorView"/> object. A new vector of the same type as this object is initialized and 
        /// returned.
        /// </summary>
        /// <param name="copyIndexingData">
        /// If true, all data of this object will be copied. If false, only the array(s) containing the values of the stored 
        /// vector entries will be copied. The new vector will reference the same indexing arrays as this one.
        /// </param>
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
        /// Calculates the Euclidian norm or 2-norm of this vector. For more see 
        /// https://en.wikipedia.org/wiki/Norm_(mathematics)#Euclidean_norm.
        /// </summary>
        double Norm2();

        /// <summary>
        /// Performs the following operation for all i: result[i] = <paramref name="scalar"/> * this[i].
        /// The resulting vector is written in a new object and then returned.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this vector.</param>
        IVector Scale(double scalar);
    }
}
