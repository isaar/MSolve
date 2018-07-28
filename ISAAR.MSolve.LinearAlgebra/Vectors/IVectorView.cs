using System;
using ISAAR.MSolve.LinearAlgebra.Reduction;

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
        IVectorView Axpy(IVectorView otherVector, double otherCoefficient);

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
        IVectorView DoEntrywise(IVectorView vector, Func<double, double, double> binaryOperation);

        /// <summary>
        /// Performs a unary operation on each entry: result[i] = <paramref name="unaryOperation"/>(this[i]).
        /// The resulting vector is written in a new object and then returned.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        IVectorView DoToAllEntries(Func<double, double> unaryOperation);
        
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
        IVectorView LinearCombination(double thisCoefficient, IVectorView otherVector, double otherCoefficient);

        /// <summary>
        /// Performs the following operation for all i: result[i] = <paramref name="scalar"/> * this[i].
        /// The resulting vector is written in a new object and then returned.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this vector.</param>
        IVectorView Scale(double scalar);
    }
}
