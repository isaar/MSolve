using System.Collections.Generic;

//TODO: provide versions of these methods that use arrays instead of dictionaries.
namespace ISAAR.MSolve.LinearAlgebra.Matrices.Builders
{
    /// <summary>
    /// Facilitates the construction of large sparse matrices, symmetric or not, usually by adding smaller submatrices. 
    /// The large matrices and their properties will be characterized as "global" in this namespace.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IGeneralMatrixBuilder: IMatrixBuilder
    {
        /// <summary>
        /// Adds the entries of <paramref name="subMatrix"/> to the corresponding entries of this matrix. The mapping between 
        /// the rows and columns of the two matrices is given by <paramref name="subRowsToGlobalRows"/> and 
        /// <paramref name="subColsToGlobalCols"/>. 
        /// </summary>
        /// <param name="subMatrix">The matrix to add to this <see cref="IGeneralMatrixBuilder"/>.</param>
        /// <param name="subRowsToGlobalRows">Mapping: row indicides of <paramref name="subMatrix"/> to row indices of 
        ///     this <see cref="IGeneralMatrixBuilder"/>.</param>
        /// <param name="subColsToGlobalCols">Mapping: column indicides of <paramref name="subMatrix"/> to column indices of 
        ///     this <see cref="IGeneralMatrixBuilder"/>.</param>
        void AddSubmatrix(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols);

        /// <summary>
        /// Adds the entries of <paramref name="subMatrix"/> to the corresponding entries of this matrix. Both matrices must be
        /// symmetric. The mapping between the rows and columns of the two matrices is given by 
        /// <paramref name="subIndicesToGlobalIndices"/> and the resulting "global" row/column indices must be symmetrically  
        /// placed around this <see cref="IGeneralMatrixBuilder"/>'s diagonal.
        /// </summary>
        /// <param name="subMatrix">The matrix to add to this <see cref="ISymmetricMatrixBuilder"/>. It must be symmetric, but 
        ///     this will not be checked explicitly.</param>
        /// <param name="subIndicesToGlobalIndices">Mapping: row/column indicides of <paramref name="subMatrix"/> to row/column 
        ///     indices of this <see cref="IGeneralMatrixBuilder"/>. The latter must be symmetrically placed around this
        ///     <see cref="IGeneralMatrixBuilder"/>'s diagonal, though this will not be checked explicitly.</param>
        void AddSubmatrixSymmetric(IIndexable2D subMatrix, IReadOnlyDictionary<int, int> subIndicesToGlobalIndices);

        /// <summary>
        /// Adds the entries of a symmetric <paramref name="subMatrix"/> to the corresponding entries of this matrix. Both 
        /// matrices must be symmetric. The mapping between the rows and columns of the two matrices is given by the pair 
        /// (<paramref name="subMatrixIndices"/>, <paramref name="globalIndices"/>) and the resulting "global" row/column   
        /// indices must be symmetrically placed around this <see cref="IGeneralMatrixBuilder"/>'s diagonal.
        /// </summary>
        /// <param name="subMatrix">
        /// The matrix to add to this <see cref="IGeneralMatrixBuilder"/>. It must be symmetric, but this will not be checked 
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
    }
}
