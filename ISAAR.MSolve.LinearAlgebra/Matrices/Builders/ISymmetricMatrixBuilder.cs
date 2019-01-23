using System.Collections.Generic;

//TODO: provide versions of these methods that use arrays instead of dictionaries.
namespace ISAAR.MSolve.LinearAlgebra.Matrices.Builders
{
    /// <summary>
    /// Facilitates the construction of large symmetric sparse matrices, usually by adding smaller submatrices.  
    /// The large matrices and their properties will be characterized as "global" in this namespace.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface ISymmetricMatrixBuilder: IMatrixBuilder
    {
        /// <summary>
        /// Adds the entries of <paramref name="subMatrix"/> to the corresponding entries of this matrix. The mapping between 
        /// the rows and columns of the two matrices is given by <paramref name="subRowsToGlobalRows"/> and 
        /// <paramref name="subColsToGlobalCols"/>. This method preprocesses each entry to find out if it is mapped to the upper 
        /// or lower triangle. As such it is safe to use, even if the user does not know which triangle the "global" row/column 
        /// indices, defined by the mapping, will end up to. If all entries are mapped to the same triangle and it is known, then
        /// <see cref="AddSubmatrixToLowerTriangle(IIndexable2D, IReadOnlyDictionary{int, int}, IReadOnlyDictionary{int, int})"/>
        /// and
        /// <see cref="AddSubmatrixToUpperTriangle(IIndexable2D, IReadOnlyDictionary{int, int}, IReadOnlyDictionary{int, int})"/>
        /// are more efficient. WARNING: Do not add (symmetric) matrices that are mapped symmetrically around this 
        /// <see cref="ISymmetricMatrixBuilder"/>'s diagonal with this method, as most entries will be added twice. Instead use
        /// <see cref="AddSubmatrixSymmetric(IIndexable2D, IReadOnlyDictionary{int, int})"/>.
        /// </summary>
        /// <param name="subMatrix">The matrix to add to this <see cref="ISymmetricMatrixBuilder"/>.</param>
        /// <param name="subRowsToGlobalRows">Mapping: row indices of <paramref name="subMatrix"/> to row indices of 
        ///     this <see cref="ISymmetricMatrixBuilder"/>.</param>
        /// <param name="subColsToGlobalCols">Mapping: column indicides of <paramref name="subMatrix"/> to column indices of 
        ///     this <see cref="ISymmetricMatrixBuilder"/>.</param>
        void AddSubmatrix(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols);

        /// <summary>
        /// Adds the entries of a symmetric <paramref name="subMatrix"/> to the corresponding entries of this matrix. The 
        /// mapping between the rows and columns of the two matrices is given by <paramref name="subIndicesToGlobalIndices"/> 
        /// and the resulting "global" row/column indices must be symmetrically placed around this 
        /// <see cref="ISymmetricMatrixBuilder"/>'s diagonal.
        /// </summary>
        /// <param name="subMatrix">The matrix to add to this <see cref="ISymmetricMatrixBuilder"/>. It must be symmetric, but 
        ///     this will not be checked explicitly.</param>
        /// <param name="subIndicesToGlobalIndices">Mapping: row/column indices of <paramref name="subMatrix"/> to row/column 
        ///     indices of this <see cref="ISymmetricMatrixBuilder"/>. The latter must be symmetrically placed around this
        ///     <see cref="ISymmetricMatrixBuilder"/>'s diagonal, though this will not be checked explicitly.</param>
        void AddSubmatrixSymmetric(IIndexable2D subMatrix, IReadOnlyDictionary<int, int> subIndicesToGlobalIndices);

        /// <summary>
        /// Adds the entries of a symmetric <paramref name="subMatrix"/> to the corresponding entries of this matrix. The 
        /// mapping between the rows and columns of the two matrices is given by the pair (<paramref name="subMatrixIndices"/>,
        /// <paramref name="globalIndices"/>) and the resulting "global" row/column indices must be symmetrically placed around  
        /// this <see cref="ISymmetricMatrixBuilder"/>'s diagonal.
        /// </summary>
        /// <param name="subMatrix">
        /// The matrix to add to this <see cref="ISymmetricMatrixBuilder"/>. It must be symmetric, but this will not be checked 
        /// explicitly.
        /// </param>
        /// <param name="subMatrixIndices">
        /// The row/column indices of the entries of <paramref name="subMatrix"/> that will be added to the global matrix. 
        /// <paramref name="subMatrixIndices"/>[i] corresponds to <paramref name="globalIndices"/>[i].
        /// </param>
        /// <param name="globalIndices">
        /// The row/column indices of the entries of the gobal matrix that will be modified. They must be symmetrically placed 
        /// around this <see cref="ISymmetricMatrixBuilder"/>'s diagonal, though this will not be checked explicitly.
        /// <paramref name="globalIndices"/>[i] corresponds to <paramref name="globalIndices"/>[i]. 
        /// </param>
        void AddSubmatrixSymmetric(IIndexable2D subMatrix, int[] subMatrixIndices, int[] globalIndices);

        /// <summary>
        /// Adds the entries of <paramref name="subMatrix"/> to the corresponding lower triangular entries of this matrix. The  
        /// mapping between the rows and columns of the two matrices is given by <paramref name="subRowsToGlobalRows"/> and 
        /// <paramref name="subColsToGlobalCols"/>. Optimized version of 
        /// <see cref="AddSubmatrix(IIndexable2D, IReadOnlyDictionary{int, int}, IReadOnlyDictionary{int, int})"/>, in case the
        /// caller is sure that all "global" row/column indices will belong to the lower triangle, though this will not be
        /// checked explicitly.
        /// </summary>
        /// <param name="subMatrix">The matrix to add to this <see cref="ISymmetricMatrixBuilder"/>.</param>
        /// <param name="subRowsToGlobalRows">Mapping: row indicides of <paramref name="subMatrix"/> to row indices of 
        ///     this <see cref="ISymmetricMatrixBuilder"/>.</param>
        /// <param name="subColsToGlobalCols">Mapping: column indices of <paramref name="subMatrix"/> to column indices of 
        ///     this <see cref="ISymmetricMatrixBuilder"/>.</param> 
        void AddSubmatrixToLowerTriangle(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols);

        /// <summary>
        /// Adds the entries of <paramref name="subMatrix"/> to the corresponding upper triangular entries of this matrix. The  
        /// mapping between the rows and columns of the two matrices is given by <paramref name="subRowsToGlobalRows"/> and 
        /// <paramref name="subColsToGlobalCols"/>. Optimized version of 
        /// <see cref="AddSubmatrix(IIndexable2D, IReadOnlyDictionary{int, int}, IReadOnlyDictionary{int, int})"/>, in case the
        /// caller is sure that all "global" row/column indices will belong to the upper triangle, though this will not be
        /// checked explicitly.
        /// </summary>
        /// <param name="subMatrix">The matrix to add to this <see cref="ISymmetricMatrixBuilder"/>.</param>
        /// <param name="subRowsToGlobalRows">Mapping: row indicides of <paramref name="subMatrix"/> to row indices of 
        ///     this <see cref="ISymmetricMatrixBuilder"/>.</param>
        /// <param name="subColsToGlobalCols">Mapping: column indices of <paramref name="subMatrix"/> to column indices of 
        ///     this <see cref="ISymmetricMatrixBuilder"/>.</param> 
        void AddSubmatrixToUpperTriangle(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols);
    }
}
