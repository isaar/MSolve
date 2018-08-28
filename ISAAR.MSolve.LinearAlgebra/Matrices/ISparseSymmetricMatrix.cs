using System.Collections.Generic;

//TODO: Perhaps the matrix itself should be the target of foreach
namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// A matrix that explicitly stores only the upper triangular non zero entries. Some zero entries may also be stored.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface ISparseSymmetricMatrix: ISymmetricMatrix, ISparseMatrix
    {
        /// <summary>
        /// Counts how many non zero entries are stored in the upper triangle of the matrix, including the diagonal. This 
        /// includes zeros that are explicitly stored.
        /// </summary>
        int CountNonZerosUpper();

        /// <summary>
        /// Iterates over the non zero entries of the upper part of the matrix, including the diagonal. This includes zeros that 
        /// are explicitly stored.
        /// </summary>
        IEnumerable<(int row, int col, double value)> EnumerateNonZerosUpper(); 
    }
}
