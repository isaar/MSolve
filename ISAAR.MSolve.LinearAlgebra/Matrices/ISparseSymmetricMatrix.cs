using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    public interface ISparseSymmetricMatrix: ISymmetricMatrix, ISparseMatrix
    {
        /// <summary>
        /// Counts the non zero entries of the upper part of the matrix, including the diagonal.
        /// </summary>
        /// <returns></returns>
        int CountNonZerosUpper();

        /// <summary>
        /// Returns the non zero entries of the upper part of the matrix, including the diagonal, as triplets (row, col, value).
        /// </summary>
        /// <returns></returns>
        IEnumerable<(int row, int col, double value)> EnumerateNonZerosSuperDiagonal(); //Perhaps the matrix itself should be the target of foreach
    }
}
