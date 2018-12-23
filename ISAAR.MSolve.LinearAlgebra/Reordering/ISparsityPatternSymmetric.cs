using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Reordering
{
    /// <summary>
    /// Describes the indices of the non-zero entries of a symmetric sparse matrix, without taking into account their values.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface ISparsityPatternSymmetric: ISparsityPattern
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
        IEnumerable<(int row, int col)> EnumerateNonZerosUpper();
    }
}
