using System;

namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    public interface IEntrywiseOperableView1D<TVectorIn, TVectorOut>
        where TVectorIn : IVectorView
        where TVectorOut : IVector
    {
        /// <summary>
        /// Performs a binary operation on each pair of entries: 
        /// result[i] = <paramref name="binaryOperation"/>(this[i], <paramref name="vector"/>[i]). 
        /// The resulting vector is written in a new object and then returned.
        /// </summary>
        /// <param name="vector">A vector with the same <see cref="IIndexable1D.Length"/> as this.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="vector"/> has different <see cref="IIndexable1D.Length"/> than this.
        /// </exception>
        TVectorOut DoEntrywise(TVectorIn vector, Func<double, double, double> binaryOperation);

        /// <summary>
        /// Performs a unary operation on each entry: result[i] = <paramref name="unaryOperation"/>(this[i]).
        /// The resulting vector is written in a new object and then returned.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        TVectorOut DoToAllEntries(Func<double, double> unaryOperation);
    }
}
