using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//TODO: integrate this interface with IMatrixBuilder and add implementations for symmetric matrices in all builders.
namespace ISAAR.MSolve.LinearAlgebra.Matrices.Builders
{
    public interface ISymmetricMatrixBuilder: IIndexable2D
    {
        /// <summary>
        /// Maps <paramref name="subMatrix"/> to "global" rows and columns of the underlying skyline matrix and then adds 
        /// the entries of <paramref name="subMatrix"/> to these global indices. The caller is respondible for the global 
        /// indices falling inside the skyline sparsity pattern. WARNING: this method adds the whole <paramref name="subMatrix"/>.
        /// If you only want to add the upper triangular part of it, use 
        /// <see cref="AddSubmatrixSymmetric(IIndexable2D, IReadOnlyDictionary{int, int})"/> instead.
        /// </summary>
        /// <param name="subMatrix"></param>
        /// <param name="subRowsToGlobalRows"></param>
        /// <param name="subColsToGlobalCols"></param>
        void AddSubmatrix(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols);

        /// <summary>
        /// Maps <paramref name="subMatrix"/> to "global" rows and columns of the underlying skyline matrix and then adds 
        /// the entries of <paramref name="subMatrix"/> to these global indices. The caller is respondible for the global 
        /// indices falling inside the sparsity pattern. Optimized version of 
        /// <see cref="AddSubmatrix(IIndexable2D, IReadOnlyDictionary{int, int}, IReadOnlyDictionary{int, int})"/>, in case the
        /// caller is sure that all global indices will be super diagonal.
        /// </summary>
        void AddSubmatrixAboveDiagonal(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols);

        /// <summary>
        /// Maps <paramref name="subMatrix"/> to "global" rows and columns of the underlying skyline matrix and then adds the
        /// entries of <paramref name="subMatrix"/> to these global indices. The caller is respondible for the global indices
        /// falling inside the sparsity pattern. Use this method if you only want to add the upper triangular part of the
        /// submatrix. Otherwise <see cref="AddSubmatrix(IIndexable2D, IReadOnlyDictionary{int, int}, IReadOnlyDictionary{int, int})"/>
        /// will add the shole submatrix, resulting in the off-diagonal entries being added twice.
        /// </summary>
        void AddSubmatrixSymmetric(IIndexable2D subMatrix, IReadOnlyDictionary<int, int> subDOFsToGlobalDOFs);

        /// <summary>
        /// Maps <paramref name="subMatrix"/> to "global" rows and columns of the underlying skyline matrix and then adds 
        /// the entries of <paramref name="subMatrix"/> to these global indices. The caller is respondible for the global 
        /// indices falling inside the sparsity pattern. Optimized version of 
        /// <see cref="AddSubmatrix(IIndexable2D, IReadOnlyDictionary{int, int}, IReadOnlyDictionary{int, int})"/>, in case the
        /// caller is sure that all global indices will be super diagonal.
        /// </summary>
        void AddSubmatrixUnderDiagonal(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols);
    }
}
