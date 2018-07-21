using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//TODO: IMatrixBuilder.AddSubmatrix() appears to be the same as ISymmetricMatrixBuilder.AddSubmatrix(). However, its sematics 
// are different. ISymmetricMatrixBuilder.AddSubmatrix() cannot be used for adding symmetric matrices, as it will add the off
// diagonal entries twice. This is due to the implementations of ISymmetricMatrixBuilder. On the other hand the implementations
// of IMatrixBuilder make it safe to pass in symmetric matrices to AddSubmatrix(), although there are more efficient methods.
namespace ISAAR.MSolve.LinearAlgebra.Matrices.Builders
{
    public interface IMatrixBuilder: IIndexable2D
    {
        new double this[int rowIdx, int colIdx] { set; }

        /// <summary>
        /// Maps <paramref name="subMatrix"/> to "global" rows and columns of the underlying skyline matrix and then adds 
        /// the entries of <paramref name="subMatrix"/> to these global indices. The caller is respondible for the global 
        /// indices falling inside the sparsity pattern. WARNING: this method adds the whole <paramref name="subMatrix"/>. 
        /// </summary>
        /// <param name="subMatrix"></param>
        /// <param name="subRowsToGlobalRows"></param>
        /// <param name="subColsToGlobalCols"></param>
        void AddSubmatrix(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols);
    }
}
