using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//TODO: perhaps the setters should be moved here too. Name them GetSubvector(), SetSubvector()
namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    /// <summary>
    /// Can return subvectors containing select entries.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface ISliceable1D: IIndexable1D
    {
        /// <summary>
        /// Returns a subvector with the entries of the original vector that correspond to the provided indices. The relative
        /// order of the entries in the returned subvector can be different than the original one. It is defined by the order 
        /// of the provided indices.
        /// </summary>
        /// <param name="indices">Which entries to return in the subvector. Constraints: 
        ///     0 &lt;= <paramref name="indices"/>[i] &lt; <see cref="IIndexable1D.Length"/>, for all i.</param>
        /// <returns></returns>
        Vector Slice(int[] indices);

        /// <summary>
        /// Returns a subvector with entries i of the original vector, such that 
        /// <paramref name="startInclusive"/> &lt;= i &lt; <paramref name="endExclusive"/> in increasing order.
        /// </summary>
        /// <param name="startInclusive">The index of the first entry that will be returned in the subvector. Constraints: 
        ///     0 &lt;= <paramref name="startInclusive"/> &lt;= <paramref name="endExclusive"/>.</param>
        /// <param name="endExclusive">The index immediately after the last entry that will be returned in the subvector.
        ///      Constraints: <paramref name="startInclusive"/> &lt;= <paramref name="endExclusive"/> &lt; 
        ///      <see cref="IIndexable1D.Length"/>.</param>
        /// <returns></returns>
        Vector Slice(int startInclusive, int endExclusive);
    }
}
