using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices.Builders
{
    public class PartialMatrixColumns
    {
        private readonly int order;
        private readonly Dictionary<int, Dictionary<int, double>> columns;

        public PartialMatrixColumns(int order)
        {
            this.order = order;
            this.columns = new Dictionary<int, Dictionary<int, double>>();
        }

        /// <summary>
        /// If the entry already exists: the new value is added to the existing one: this[rowIdx, colIdx] += value. 
        /// Otherwise the entry is inserted with the new value: this[rowIdx, colIdx] = value.
        /// Use this method instead of this[rowIdx, colIdx] += value, as it is an optimized version.
        /// Since both sub-diagonal and super-diagonal entries of the symmetric matrix are stored separately, it is safe to
        /// process both (<paramref name="rowIdx"/>, <paramref name="colIdx"/>) and 
        /// (<paramref name="colIdx"/>, <paramref name="rowIdx"/>).
        /// </summary>
        /// <param name="rowIdx"></param>
        /// <param name="colIdx"></param>
        /// <param name="value"></param>
        public void AddToEntry(int rowIdx, int colIdx, double value)
        {
            if (columns.TryGetValue(colIdx, out Dictionary<int, double> wholeColumn)) // The column exists. Mutate it.
            {
                // If the entry did not exist, oldValue = default(double) = 0.0 and adding it should not reduce precision.
                wholeColumn.TryGetValue(rowIdx, out double oldValue);
                wholeColumn[rowIdx] = value + oldValue;
                //The Dictionary wholeColumn is indexed twice in both cases. Is it possible to only index it once?
            }
            else // The column did not exist. Create a new one with the desired entry.
            {
                wholeColumn = new Dictionary<int, double>();
                wholeColumn.Add(rowIdx, value);
                columns.Add(colIdx, wholeColumn);
            }
        }

        /// <summary>
        /// Returns the whole column with index = <paramref name="colIdx"/>. However, the entries with row index that belongs in
        /// <paramref name="tabooRows"/> will be set to 0. Actually they will not be included in the sparsity pattern of the
        /// returned <see cref="SparseVector"/>. Note that the length of the returned vector is equal to <see cref="NumRows"/>. 
        /// Since the matrix is symmetric, row = column. 
        /// </summary>
        /// <param name="colIdx">The index of the column to return.</param>
        /// <param name="tabooRows">The entries of the returned column vector at the indices <paramref name="tabooRows"/> will 
        ///     be equal to 0.</param>
        /// <returns></returns>
        public SparseVector GetColumnWithoutRows(int colIdx, HashSet<int> tabooRows)
        {
            bool exists = columns.TryGetValue(colIdx, out Dictionary<int, double> wholeColumn); //TODO: Should I just let the Dictionary indexer throw?
            if (!exists) throw new IndexOutOfRangeException($"Column {colIdx} is not stored.");
            var result = new SortedDictionary<int, double>();
            foreach (var rowVal in columns[colIdx])
            {
                if (!tabooRows.Contains(rowVal.Key)) result.Add(rowVal.Key, rowVal.Value);
            }
            return SparseVector.CreateFromDictionary(order, result);
        }

        public SparseVector GetColumnWithSelectedRows(int colIdx, HashSet<int> wantedRows)
        {
            //TODO: perhaps I should pass a sorted set, iterate it and create the raw arrays directly. Or sort the passed set.
            bool exists = columns.TryGetValue(colIdx, out Dictionary<int, double> wholeColumn); //TODO: Should I just let the Dictionary indexer throw?
            if (!exists) throw new IndexOutOfRangeException($"Column {colIdx} is not stored.");
            var result = new SortedDictionary<int, double>();
            foreach (var rowVal in columns[colIdx])
            {
                if (wantedRows.Contains(rowVal.Key)) result.Add(rowVal.Key, rowVal.Value);
            }
            return SparseVector.CreateFromDictionary(order, result);
        }

        //TODO: this is probably not needed as the caller could approach identity columns more efficiently
        public void SetColumnToIdentity(int colIdx)
        {
            var identityCol = new Dictionary<int, double>();
            identityCol[colIdx] = 1.0;
            columns[colIdx] = identityCol;
        }
    }
}
