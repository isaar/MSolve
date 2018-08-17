using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;

//TODO: Perhaps the matrix itself should be the target of foreach. 
namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// A matrix that does not explicitly store all/most of its zero entries.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface ISparseMatrix: IIndexable2D
    {
        /// <summary>
        /// Counts how many non zero entries are stored in the matrix. This includes zeros that are explicitly stored.
        /// </summary>
        int CountNonZeros();

        /// <summary>
        /// Iterates over the non zero entries of the matrix. This includes zeros that are explicitly stored.
        /// </summary>
        IEnumerable<(int row, int col, double value)> EnumerateNonZeros();

        /// <summary>
        /// Returns a data transfer object that contains the internal data structures of this matrix.
        /// </summary>
        SparseFormat GetSparseFormat();
    }
}
