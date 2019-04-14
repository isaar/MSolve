using System;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    public interface IEntrywiseOperable2D<TMatrixIn>
        where TMatrixIn : IMatrixView
    {
        /// <summary>
        /// Performs a binary operation on each pair of entries:  
        /// this[i, j] = <paramref name="binaryOperation"/>(this[i,j], <paramref name="matrix"/>[i,j]).
        /// The resulting matrix overwrites the entries of this.
        /// </summary>
        /// <param name="matrix">
        /// A matrix with the same <see cref="IIndexable2D.NumRows"/> and <see cref="IIndexable2D.NumColumns"/> as this.
        /// </param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="matrix"/> has different <see cref="IIndexable2D.NumRows"/> or 
        /// <see cref="IIndexable2D.NumColumns"/> than this.
        /// </exception>
        /// <exception cref="Exceptions.PatternModifiedException">
        /// Thrown if an entry this[i, j] needs to be overwritten, but that is not permitted by the matrix storage format.
        /// </exception>
        void DoEntrywiseIntoThis(TMatrixIn matrix, Func<double, double, double> binaryOperation);

        /// <summary>
        /// Performs a unary operation on each entry: this[i] = <paramref name="unaryOperation"/>(this[i, j]).
        /// he resulting matrix overwrites the entries of this.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        /// <exception cref="Exceptions.PatternModifiedException">
        /// Thrown if an entry this[i, j] needs to be overwritten, but that is not permitted by the matrix storage format.
        /// </exception>
        void DoToAllEntriesIntoThis(Func<double, double> unaryOperation);
    }
}
