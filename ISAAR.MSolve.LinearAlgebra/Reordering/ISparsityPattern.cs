using System.Collections.Generic;

namespace ISAAR.MSolve.LinearAlgebra.Reordering
{
    /// <summary>
    /// Describes the indices of the non-zero entries of a sparse matrix, without taking into account their values.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface ISparsityPattern
    {
        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        int NumColumns { get; }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        int NumRows { get; }

        /// <summary>
        /// Counts how many non zero entries are stored in the matrix. This includes zeros that are explicitly stored.
        /// </summary>
        int CountNonZeros();

        /// <summary>
        /// Iterates over the non zero entries of the matrix. This includes zeros that are explicitly stored.
        /// </summary>
        IEnumerable<(int row, int col)> EnumerateNonZeros();

        /// <summary>
        /// Returns true if the matrix entry at (<paramref name="rowIdx"/>, <paramref name="colIdx"/>) is non-zero or it is 
        /// equal to 0, but explicitly stored (a.k.a. non-structural zero). Otherwise returns false.
        /// </summary>
        /// <param name="rowIdx">The row index: 0 &lt;= <paramref name="rowIdx"/> &lt; <see cref="NumRows"/>.</param>
        /// <param name="colIdx">The column index: 0 &lt;= <paramref name="colIdx"/> &lt; <see cref="NumColumns"/>.</param>
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="rowIdx"/> or <paramref name="colIdx"/> violate 
        ///     the described constraints.</exception>
        bool IsNonZero(int rowIdx, int colIdx);
    }
}
