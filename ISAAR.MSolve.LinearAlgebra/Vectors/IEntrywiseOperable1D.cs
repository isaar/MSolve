using System;

//TODO: Should Axpy() and LinearCombination be here as well? 
//TODO: Should this extend IEntrywiseOperableView1D?
namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    public interface IEntrywiseOperable1D<TVectorIn>
        where TVectorIn : IVectorView
    {
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
        void DoEntrywiseIntoThis(TVectorIn otherVector, Func<double, double, double> binaryOperation);

        /// <summary>
        /// Performs a unary operation on each entry: this[i] = <paramref name="unaryOperation"/>(this[i]).
        /// The resulting vector overwrites the entries of this.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        /// <exception cref="Exceptions.PatternModifiedException">
        /// Thrown if an entry this[i] needs to be overwritten, but that is not permitted by the vector storage format.
        /// </exception> 
        void DoToAllEntriesIntoThis(Func<double, double> unaryOperation);
    }
}
