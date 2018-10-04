using System;
using System.Collections.Generic;

namespace ISAAR.MSolve.Numerical.Commons
{
    /// <summary>
    /// Interface for table data structures, which associate ordered pairs (row, column) with values. Contrary to matrices, not 
    /// all (row, column) pairs of a table need to be associated with values.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    /// <typeparam name="TRow"></typeparam>
    /// <typeparam name="TColumn"></typeparam>
    /// <typeparam name="TValue"></typeparam>
    public interface ITable<TRow, TColumn, TValue> : IEnumerable<(TRow row, TColumn col, TValue val)>
    {
        /// <summary>
        /// The total number of entries in this table. 
        /// </summary>
        int EntryCount { get; }

        /// <summary>
        /// The value associated with the pair (<paramref name="row"/>, <paramref name="col"/>).
        /// </summary>
        /// <param name="row"></param>
        /// <param name="col"></param>
        /// <exception cref="KeyNotFoundException">Thrown if the pair (<paramref name="row"/>, <paramref name="col"/>) has not 
        ///     been associated with a value.</exception>
        TValue this[TRow row, TColumn col] { get; }

        /// <summary>
        /// True if the pair (<paramref name="row"/>, <paramref name="col"/>) has already been associated with a value.
        /// </summary>
        /// <param name="row"></param>
        /// <param name="col"></param>
        bool Contains(TRow row, TColumn col);

        /// <summary>
        /// Returns a collection of the <see cref="TColumn"/> indices of all (<paramref name="row"/>, <see cref="TColumn"/>) 
        /// pairs in this object.
        /// </summary>
        /// <param name="row"></param>
        IEnumerable<TColumn> GetColumnsOfRow(TRow row);

        /// <summary>
        /// Returns a collection of all the <see cref="TRow"/> indices in this object.
        /// </summary>
        IEnumerable<TRow> GetRows();

        /// <summary>
        /// Returns a collection of all <see cref="TValue"/> values associated with the (<paramref name="row"/>, 
        /// <see cref="TColumn"/>) pairs in this object.
        /// </summary>
        /// <param name="row"></param>
        IEnumerable<TValue> GetValuesOfRow(TRow row);

        /// <summary>
        /// Returns true if the pair (<paramref name="row"/>, <paramref name="col"/>) has been associated with a value, which is  
        /// also returned as an output parameter <paramref name="value"/>. Otherwise returns false instead of throwing any 
        /// exception.
        /// </summary>
        /// <param name="row"></param>
        /// <param name="col"></param>
        /// <param name="value"></param>
        bool TryGetValue(TRow row, TColumn col, out TValue value);
    }
}