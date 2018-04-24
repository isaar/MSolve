using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.LinearAlgebra.Matrices.Builders
{
    public interface IMatrixBuilder: IIndexable2D
    {
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
