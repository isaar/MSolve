using System;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Perhaps Addition, Subtraction and Scaling must be done without using delegates, for performance
//TODO: perhaps I should return IMatrixView instead of IMatrix. By returning IMatrixView I can have classes that only implement
//      IMatrixView. On the other hand, I cannot mutate the returned type, so its usefulness is limited.
namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// It supports common operations that do not mutate the underlying matrix. If you need to store a matrix and then pass it
    /// around or allow acceess to it, consider using this interface instead of <see cref="Matrix"/> for extra safety.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IMatrixView: IIndexable2D, IReducible
    {
        /// <summary>
        /// Performs the following operation for all (i, j):
        /// result[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
        /// Optimized version of <see cref="IMatrixView.DoEntrywise(IMatrixView, Func{double, double, double})"/> and 
        /// <see cref="IMatrixView.LinearCombination(double, IMatrixView, double)"/>. Named after BLAS axpy (y = a * x plus y).
        /// The resulting matrix is written in a new object and then returned.
        /// </summary>
        /// <param name="otherMatrix">A matrix with the same <see cref="IIndexable2D.NumRows"/> and 
        ///     <see cref="IIndexable2D.NumColumns"/> as this.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     <see cref="IIndexable2D.NumRows"/> or <see cref="IIndexable2D.NumColumns"/> than this.</exception>
        IMatrix Axpy(IMatrixView otherMatrix, double otherCoefficient);


        /// <summary>
        /// Performs a binary operation on each pair of entries: 
        /// result[i, j] = <paramref name="binaryOperation"/>(this[i, j], <paramref name="matrix"/>[i]). 
        /// The resulting matrix is written in a new object and then returned.
        /// </summary>
        /// <param name="matrix">A matrix with the same <see cref="IIndexable2D.NumRows"/> and 
        ///     <see cref="IIndexable2D.NumColumns"/> as this.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if <paramref name="matrix"/> has different 
        ///     <see cref="IIndexable2D.NumRows"/> or <see cref="IIndexable2D.NumColumns"/> than this.</exception>
        IMatrix DoEntrywise(IMatrixView matrix, Func<double, double, double> binaryOperation);

        /// <summary>
        /// Performs a unary operation on each entry: result[i] = <paramref name="unaryOperation"/>(this[i, j]).
        /// The resulting matrix is written in a new object and then returned.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        IMatrix DoToAllEntries(Func<double, double> unaryOperation);


        /// <summary>
        /// Performs the following operation for all (i, j):
        /// result[i, j] = <paramref name="thisCoefficient"/> * this[i, j] + <paramref name="otherCoefficient"/> * 
        /// <paramref name="otherMatrix"/>[i, j]. 
        /// Optimized version of <see cref="IMatrixView.DoEntrywise(IMatrixView, Func{double, double, double})"/>. 
        /// The resulting matrix is written in a new object and then returned.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this.</param>
        /// <param name="otherMatrix">A matrix with the same <see cref="IIndexable2D.NumRows"/> and 
        ///     <see cref="IIndexable2D.NumColumns"/> as this.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     <see cref="IIndexable2D.NumRows"/> or <see cref="IIndexable2D.NumColumns"/> than this.</exception>
        IMatrix LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient);

        /// <summary>
        /// Performs the matrix-matrix multiplication: oper(<paramref name="other"/>) * oper(this).
        /// </summary>
        /// <param name="other">A matrix such that the <see cref="IIndexable2D.NumColumns"/> of oper(<paramref name="other"/>) 
        ///     are equal to the <see cref="IIndexable2D.NumRows"/> of oper(this).</param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
        /// <param name="transposeOther">If true, oper(<paramref name="other"/>) = transpose(<paramref name="other"/>). 
        ///     Otherwise oper(<paramref name="other"/>) = <paramref name="other"/>.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if oper(<paramref name="otherMatrix"/>) has 
        ///     different <see cref="IIndexable2D.NumColumns"/> than the <see cref="IIndexable2D.NumRows"/> of 
        ///     oper(this).</exception>
        Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false);

        /// <summary>
        /// Performs the matrix-matrix multiplication: oper(this) * oper(<paramref name="other"/>).
        /// </summary>
        /// <param name="other">A matrix such that the <see cref="IIndexable2D.NumRows"/> of oper(<paramref name="other"/>) 
        ///     are equal to the <see cref="IIndexable2D.NumColumns"/> of oper(this).</param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
        /// <param name="transposeOther">If true, oper(<paramref name="other"/>) = transpose(<paramref name="other"/>). 
        ///     Otherwise oper(<paramref name="other"/>) = <paramref name="other"/>.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if oper(<paramref name="otherMatrix"/>) has 
        ///     different <see cref="IIndexable2D.NumRows"/> than the <see cref="IIndexable2D.NumColumns"/> of 
        ///     oper(this).</exception>
        Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false);

        /// <summary>
        /// Performs the matrix-vector multiplication: oper(this) * <paramref name="vector"/>.
        /// To multiply this * columnVector, set <paramref name="transposeThis"/> to false.
        /// To multiply rowVector * this, set <paramref name="transposeThis"/> to true.
        /// The resulting vector will be written in a new vector and returned.
        /// </summary>
        /// <param name="vector">A vector with <see cref="IIndexable1D.Length"/> being equal to the 
        ///     <see cref="IIndexable2D.NumColumns"/> of oper(this).</param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if the <see cref="IIndexable1D.Length"/> of
        ///     <paramref name="vector"/> is different than the <see cref="IIndexable2D.NumColumns"/> of oper(this).</exception>
        IVector Multiply(IVectorView vector, bool transposeThis = false);

        /// <summary>
        /// Performs the matrix-vector multiplication: <paramref name="rhsVector"/> = oper(this) * <paramref name="vector"/>.
        /// To multiply this * columnVector, set <paramref name="transposeThis"/> to false.
        /// To multiply rowVector * this, set <paramref name="transposeThis"/> to true.
        /// The resulting vector will overwrite the entries of <paramref name="rhsVector"/>.
        /// </summary>
        /// <param name="lhsVector">
        /// The vector that will be multiplied by this matrix. It sits on the left hand side of the equation y = oper(A) * x.
        /// Constraints: <paramref name="lhsVector"/>.<see cref="IIndexable1D.Length"/> 
        /// == oper(this).<see cref="IIndexable2D.NumColumns"/>.
        /// </param>
        /// <param name="rhsVector">
        /// The vector that will be overwritten by the result of the multiplication. It sits on the right hand side of the 
        /// equation y = oper(A) * x. Constraints: <paramref name="lhsVector"/>.<see cref="IIndexable1D.Length"/> 
        /// == oper(this).<see cref="IIndexable2D.NumRows"/>.
        /// </param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if the <see cref="IIndexable1D.Length"/> of <paramref name="lhsVector"/> or <paramref name="rhsVector"/> 
        /// violate the described contraints.
        /// </exception>
        /// <exception cref="Exceptions.PatternModifiedException">
        /// Thrown if the storage format of <paramref name="rhsVector"/> does not support overwritting the entries that this 
        /// method will try to.
        /// </exception>
        void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis);

        /// <summary>
        /// Performs the following operation for all (i, j): result[i, j] = <paramref name="scalar"/> * this[i, j].
        /// The resulting matrix is written in a new object and then returned.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this matrix.</param>
        IMatrix Scale(double scalar);

        /// <summary>
        /// Returns a matrix that is transpose to this: result[i, j] = this[j, i]. The entries will be explicitly copied. Some
        /// implementations of <see cref="IMatrixView"/> may offer more efficient transpositions, that do not copy the entries.
        /// If the transposed matrix will be used only for multiplications, <see cref="MultiplyLeft(IMatrixView, bool, bool)"/>,
        /// <see cref="MultiplyRight(IMatrixView, bool, bool)"/> and <see cref="Multiply(IVectorView, bool)"/> are more 
        /// effient generally.
        /// </summary>
        IMatrix Transpose(); //TODO: perhaps this should default to not copying the entries, if possible.
    }
}
