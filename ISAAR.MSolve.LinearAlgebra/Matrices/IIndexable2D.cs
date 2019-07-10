namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// A matrix that supports indexing These are the most basic operations all matrix classes must
    /// implement. As such, it can be used for matrix formats that do not support linear algebra operations, such as DOKs and 
    /// other builders.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IIndexable2D : IBounded2D
    {
        /// <summary>
        /// The entry with row index = rowIdx and column index = colIdx. 
        /// </summary>
        /// <param name="rowIdx">The row index: 0 &lt;= <paramref name="rowIdx"/> &lt; <see cref="NumRows"/>.</param>
        /// <param name="colIdx">The column index: 0 &lt;= <paramref name="colIdx"/> &lt; <see cref="NumColumns"/>.</param>
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="rowIdx"/> or <paramref name="colIdx"/> violate 
        ///     the described constraints.</exception>
        double this[int rowIdx, int colIdx] { get; }

        /// <summary>
        /// Returns true if 1) this and <paramref name="other"/> have the same <see cref="NumRows"/> and <see cref="NumColumns"/>
        /// 2) this[i, j] - <paramref name="other"/>[i, j] is within the acceptable <paramref name="tolerance"/> for all (i, j). 
        /// </summary>
        /// <param name="other">The other matrix that this matrix will be compared to.</param>
        /// <param name="tolerance">The entries at (i, j) of the two matrices will be considered equal, if
        ///     (<paramref name="other"/>[i, j] - this[i, j]) / this[i, j] &lt;= <paramref name="tolerance"/>. Setting 
        ///     <paramref name="tolerance"/> = 0, will check if these entries are exactly the same.</param>
        bool Equals(IIndexable2D other, double tolerance = 1e-13);
    }
}
