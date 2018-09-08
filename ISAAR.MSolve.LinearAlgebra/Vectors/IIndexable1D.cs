namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    /// <summary>
    /// A vector that supports indexing and dimension querying. These are the most basic operations all vector classes must
    /// implement.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IIndexable1D
    {
        /// <summary>
        /// The number of entries of the vector. 
        /// </summary>
        int Length { get; }

        /// <summary>
        /// Returns the entry at <paramref name="index"/>. 
        /// </summary>
        /// <param name="index">The index of the entry to return. Constraints: 
        ///     0 &lt;= <paramref name="index"/> &lt; <see cref="Length"/>.</param>
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="index"/> violates the described constraint.
        ///     </exception>
        double this[int index] { get; }

        /// <summary>
        /// Returns true if 1) this and <paramref name="other"/> have the same <see cref="Length"/> and 
        /// 2) this[i] - <paramref name="other"/>[i] is within the acceptable <paramref name="tolerance"/> for all i. 
        /// </summary>
        /// <param name="other">The other vector that this vector will be compared to.</param>
        /// <param name="tolerance">The entries at index i of the two vectors will be considered equal, if
        ///     (<paramref name="other"/>[i] - this[i]) / this[i] &lt;= <paramref name="tolerance"/>. Setting 
        ///     <paramref name="tolerance"/> = 0, will check if these entries are exactly the same.</param>
        bool Equals(IIndexable1D other, double tolerance = 1e-13);
    }
}
