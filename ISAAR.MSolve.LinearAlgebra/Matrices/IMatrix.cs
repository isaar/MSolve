using System;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Operations specified by this interface modify the matrix. Therefore it is possible that they may throw exceptions if they 
    /// are used on sparse or triangular storage matrix formats.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IMatrix: IMatrixView
    {
        /// <summary>
        /// Performs the following operation for all (i, j):
        /// this[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
        /// Optimized version of <see cref="DoEntrywiseIntoThis(IMatrixView, Func{double, double, double})"/> and 
        /// <see cref="LinearCombinationIntoThis(double, IMatrixView, double)"/>. Named after BLAS axpy (y = a*x plus y). 
        /// The resulting matrix overwrites the entries of this.
        /// </summary>
        /// <param name="other">A matrix with the same <see cref="IIndexable2D.NumRows"/> and 
        ///     <see cref="IIndexable2D.NumColumns"/> as this.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     <see cref="IIndexable2D.NumRows"/> or <see cref="IIndexable2D.NumColumns"/> than this.</exception>
        /// <exception cref="Exceptions.PatternModifiedException">Thrown if an entry this[i, j] needs to be overwritten, but that 
        ///     is not permitted by the matrix storage format.</exception>
        void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient);

        /// <summary>
        /// Sets all entries to 0. For sparse or block matrices: the indexing arrays will not be mutated. Therefore the sparsity  
        /// pattern will be preserved. The non-zero entries will be set to 0, but they will still be stored explicitly. 
        /// </summary>
        void Clear();

        /// <summary>
        /// Performs a binary operation on each pair of entries:  
        /// this[i, j] = <paramref name="binaryOperation"/>(this[i,j], <paramref name="matrix"/>[i,j]).
        /// The resulting matrix overwrites the entries of this.
        /// </summary>
        /// <param name="matrix">A matrix with the same <see cref="IIndexable2D.NumRows"/> and 
        ///     <see cref="IIndexable2D.NumColumns"/> as this.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if <paramref name="matrix"/> has different 
        ///     <see cref="IIndexable2D.NumRows"/> or <see cref="IIndexable2D.NumColumns"/> than this.</exception>
        /// <exception cref="Exceptions.PatternModifiedException">Thrown if an entry this[i, j] needs to be overwritten, but that 
        ///     is not permitted by the matrix storage format.</exception>
        void DoEntrywiseIntoThis(IMatrixView matrix, Func<double, double, double> binaryOperation);

        /// <summary>
        /// Performs a unary operation on each entry: this[i] = <paramref name="unaryOperation"/>(this[i, j]).
        /// he resulting matrix overwrites the entries of this.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        /// <exception cref="Exceptions.PatternModifiedException">Thrown if an entry this[i, j] needs to be overwritten, but that 
        ///     is not permitted by the matrix storage format.</exception>
        void DoToAllEntriesIntoThis(Func<double, double> unaryOperation);

        /// <summary>
        /// Performs the following operation for all (i, j):
        /// this[i, j] = <paramref name="thisCoefficient"/> * this[i, j] + <paramref name="otherCoefficient"/> * 
        /// <paramref name="otherMatrix"/>[i, j]. 
        /// Optimized version of <see cref="DoEntrywiseIntoThis(IMatrixView, Func{double, double, double})"/>.
        /// The resulting matrix overwrites the entries of this.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this.</param>
        /// <param name="otherMatrix">A matrix with the same <see cref="IIndexable2D.NumRows"/> and 
        ///     <see cref="IIndexable2D.NumColumns"/> as this.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     <see cref="IIndexable2D.NumRows"/> or <see cref="IIndexable2D.NumColumns"/> than this.</exception>
        /// <exception cref="Exceptions.PatternModifiedException">Thrown if an entry this[i, j] needs to be overwritten, but that 
        ///     is not permitted by the matrix storage format.</exception>
        void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient);

        /// <summary>
        /// Performs the following operation for all (i, j): this[i, j] = <paramref name="scalar"/> * this[i, j].
        /// The resulting matrix overwrites the entries of this.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this matrix.</param>
        void ScaleIntoThis(double scalar);

        /// <summary>
        /// Setter that will work as expected for general dense matrices. For sparse matrices it will throw a 
        /// <see cref="Exceptions.SparsityPatternModifiedException"/> if a structural zero entry is written to.
        /// For symmetric matrices, this will set both (<paramref name="rowIdx"/>, <paramref name="colIdx"/>) and 
        /// (<paramref name="colIdx"/>, <paramref name="rowIdx"/>).
        /// </summary>
        /// <param name="rowIdx">The row index: 0 &lt;= rowIdx &lt; <see cref="IIndexable2D.NumRows"/></param>
        /// <param name="colIdx">The column index: 0 &lt;= colIdx &lt; <see cref="IIndexable2D.NumColumns"/></param>
        /// <param name="value">The new value of this[<paramref name="rowIdx"/>, <paramref name="colIdx"/>].</param>
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="rowIdx"/> or <paramref name="colIdx"/> violate 
        ///     the described constraints.</exception>
        /// <exception cref="Exceptions.SparsityPatternModifiedException"> Thrown if a structural zero entry of a sparse matrix 
        ///     format is written to.</exception>
        void SetEntryRespectingPattern(int rowIdx, int colIdx, double value);
    }
}
