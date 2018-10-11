using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Numerical.Commons
{
    /// <summary>
    /// Interface for table data structures, which associate ordered triplets (key1, key2, key3) with values. Not all triplets
    /// of a table need to be associated with values.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    /// <typeparam name="TKey1"></typeparam>
    /// <typeparam name="TKey2"></typeparam>
    /// <typeparam name="TKey3"></typeparam>
    /// <typeparam name="TValue"></typeparam>
    public interface ITable3D<TKey1, TKey2, TKey3, TValue> : IEnumerable<(TKey1 key1, TKey2 key2, TKey3 key3, TValue val)>
    {
        /// <summary>
        /// The total number of entries in this table. 
        /// </summary>
        int EntryCount { get; }

        /// <summary>
        /// The value associated with the triplet (<paramref name="key1"/>, <paramref name="key2"/>, <paramref name="key3"/>).
        /// </summary>
        /// <param name="key1"></param>
        /// <param name="key2"></param>
        /// <param name="key3"></param>
        /// <exception cref="KeyNotFoundException">Thrown if the triplet 
        ///     (<paramref name="row"/>, <paramref name="key2"/>, <paramref name="key3"/>) has not been associated with a
        ///     value.</exception>
        TValue this[TKey1 key1, TKey2 key2, TKey3 key3] { get; }

        /// <summary>
        /// True if the triplet (<paramref name="key1"/>, <paramref name="key2"/>, <paramref name="key3"/>) has already been 
        /// associated with a value.
        /// </summary>
        /// <param name="key1"></param>
        /// <param name="key2"></param>
        /// <param name="key3"></param>
        bool Contains(TKey1 key1, TKey2 key2, TKey3 key3);

        /// <summary>
        /// Returns true if the triplet (<paramref name="key1"/>, <paramref name="key2"/>, <paramref name="key3"/>) has been 
        /// associated with a value, which is also returned as an output parameter <paramref name="value"/>. Otherwise returns 
        /// false instead of throwing any exception.
        /// </summary>
        /// <param name="key1"></param>
        /// <param name="key2"></param>
        /// <param name="key3"></param>
        /// <param name="value"></param>
        bool TryGetValue(TKey1 key1, TKey2 key2, TKey3 key3, out TValue value);
    }
}
