using System;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    public interface IEntrywiseOperableView2D<TMatrixIn, TMatrixOut>
        where TMatrixIn : IMatrixView
        where TMatrixOut : IMatrix
    {
        /// <summary>
        /// Performs a binary operation on each pair of entries: 
        /// result[i, j] = <paramref name="binaryOperation"/>(this[i, j], <paramref name="matrix"/>[i]). 
        /// The resulting matrix is written in a new object and then returned.
        /// </summary>
        /// <param name="matrix">
        /// A matrix with the same <see cref="IIndexable2D.NumRows"/> and <see cref="IIndexable2D.NumColumns"/> as this.
        /// </param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="matrix"/> has different <see cref="IIndexable2D.NumRows"/> or 
        /// <see cref="IIndexable2D.NumColumns"/> than this.
        /// </exception>
        TMatrixOut DoEntrywise(TMatrixIn matrix, Func<double, double, double> binaryOperation);

        /// <summary>
        /// Performs a unary operation on each entry: result[i] = <paramref name="unaryOperation"/>(this[i, j]).
        /// The resulting matrix is written in a new object and then returned.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        TMatrixOut DoToAllEntries(Func<double, double> unaryOperation);
    }
}
