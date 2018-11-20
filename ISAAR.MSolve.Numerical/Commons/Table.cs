using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;

//TODO: Add slicing features
//TODO: Perhaps a Table.Builder would be better, since EntryCount can be O(1)
//TODO: If we know that the rows are ALL the nodes of the model and they are numberd contiguously, we can implement the rows 
//      as an array instead of Dictionary, in order to speed up things.
//TODO: Many of the methods will be called multiple times with the same row. Provide accessors for a TableRow object to let the
//      client avoid looking up the row all the time.
//TODO: cache the last row looked up (and its Dictionary<TColumn, TValue>) to improve performance for consecutive look ups 
namespace ISAAR.MSolve.Numerical.Commons
{
    /// <summary>
    /// Basic table data structure, which associates ordered pairs (row, column) with values. Contrary to matrices, not 
    /// all (row, column) pairs of a table need to be associated with values.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    /// <typeparam name="TRow"></typeparam>
    /// <typeparam name="TColumn"></typeparam>
    /// <typeparam name="TValue"></typeparam>
    public class Table<TRow, TColumn, TValue> : ITable<TRow, TColumn, TValue>
    {
        protected const int defaultInitialCapacity = 1; //There will be at least 1. TODO: perhaps get this from Dictionary class.
        protected readonly int initialCapacityForEachDim;
        protected readonly Dictionary<TRow, Dictionary<TColumn, TValue>> data;

        public Table(int initialCapacityForEachDim = defaultInitialCapacity)
        {
            this.initialCapacityForEachDim = initialCapacityForEachDim;
            this.data = new Dictionary<TRow, Dictionary<TColumn, TValue>>(initialCapacityForEachDim);
        }

        protected Table(Dictionary<TRow, Dictionary<TColumn, TValue>> data,
            int initialCapacityForEachDim = defaultInitialCapacity)
        {
            this.initialCapacityForEachDim = initialCapacityForEachDim;
            this.data = data;
        }

        public int EntryCount //TODO: perhaps this should be cached somehow
        {
            get
            {
                int count = 0;
                foreach (var wholeRow in data) count += wholeRow.Value.Count;
                return count;
            }
        }

        public TValue this[TRow row, TColumn col]
        {
            get => data[row][col];

            set
            {
                bool containsRow = data.TryGetValue(row, out Dictionary<TColumn, TValue> wholeRow);
                if (!containsRow)
                {
                    wholeRow = new Dictionary<TColumn, TValue>(initialCapacityForEachDim);
                    data.Add(row, wholeRow);
                }
                wholeRow[col] = value; // This allows changing the value after an entry has been added.
                // The code below was used to prevent changes after an entry has been added, but that makes reordering difficult.
                //if (wholeRow.ContainsKey(col))
                //{
                //    throw new ArgumentException("The entry (row, column) = (" + row + ", "
                //        + col + ") already exists.");
                //}
                //else wholeRow.Add(col, value);
            }
        }

        public void Clear() => data.Clear();

        public bool Contains(TRow row, TColumn col)
        {
            bool containsRow = data.TryGetValue(row, out Dictionary<TColumn, TValue> wholeRow);
            if (!containsRow) return false;
            else return wholeRow.ContainsKey(col);
        }

        public IEnumerator<(TRow row, TColumn col, TValue val)> GetEnumerator()
        {
            foreach (var wholeRow in data)
            {
                foreach (var colValPair in wholeRow.Value)
                {
                    yield return (wholeRow.Key, colValPair.Key, colValPair.Value);
                }
            }
        }

        IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

        public IEnumerable<TColumn> GetColumnsOfRow(TRow row)
        {
            bool containsRow = data.TryGetValue(row, out Dictionary<TColumn, TValue> wholeRow);
            if (containsRow) return wholeRow.Keys;
            else return Enumerable.Empty<TColumn>();
        }

        public IEnumerable<TRow> GetRows() => data.Keys;

        public IEnumerable<TValue> GetValuesOfRow(TRow row)
        {
            bool containsRow = data.TryGetValue(row, out Dictionary<TColumn, TValue> wholeRow);
            if (containsRow) return wholeRow.Values;
            else return Enumerable.Empty<TValue>();
        }

        public void ModifyValues(Func<TValue, TValue> unaryOperation)
        {
            //TODO: perhaps I should create a new table and replace the existing one once finished.

            foreach (Dictionary<TColumn, TValue> rowData in data.Values)
            {
                // Create a temporary collection for iterating, while modifying the actual one
                var rowDataAsList = new List<KeyValuePair<TColumn, TValue>>(rowData);
                
                foreach (var colValPair in rowDataAsList)
                {
                    rowData[colValPair.Key] = unaryOperation(colValPair.Value);
                }
            }
        }

        public override string ToString()
        {
            StringBuilder builder = new StringBuilder();
            foreach (var wholeRow in this)
            {
                builder.Append(wholeRow.Item1);
                builder.Append(" , ");
                builder.Append(wholeRow.Item2);
                builder.Append(" , ");
                builder.Append(wholeRow.Item3);
                builder.Append('\n');
            }
            return builder.ToString();
        }

        /// <summary>
        /// Adds the specified (<paramref name="row"/>, <paramref name="col"/>, <paramref name="value"/>) entry to the table. 
        /// Returns true, if the insertion was successful, or false, if the table already contained the specified entry.
        /// </summary>
        /// <param name="row">The </param>
        /// <param name="col"></param>
        /// <param name="value"></param>
        public bool TryAdd(TRow row, TColumn col, TValue value)
        {
            bool containsRow = data.TryGetValue(row, out Dictionary<TColumn, TValue> wholeRow);
            if (containsRow)
            {
                //TODO: the following try clause can be replaced by the more efficient "return wholeRow.TryAdd(col, value)", once
                // it is available in .NET Standard. Unfortunately I cannot reference a .Net Core 2.1 project from .Net Standard
                try
                {
                    wholeRow.Add(col, value);
                    return true;
                }
                catch (ArgumentException)
                {
                    return false;
                }
            }
            else
            {
                wholeRow = new Dictionary<TColumn, TValue>(initialCapacityForEachDim);
                data.Add(row, wholeRow);
                wholeRow.Add(col, value);
                return true;
            }
        }

        public bool TryGetDataOfRow(TRow row, out IReadOnlyDictionary<TColumn, TValue> columnValuePairs)
        {
            bool exists = data.TryGetValue(row, out Dictionary<TColumn, TValue> rowData);
            columnValuePairs = rowData;
            return exists;
        }

        public bool TryGetValue(TRow row, TColumn col, out TValue value)
        {
            bool containsRow = data.TryGetValue(row, out Dictionary<TColumn, TValue> wholeRow);
            if (!containsRow)
            {
                value = default(TValue);
                return false;
            }
            else return wholeRow.TryGetValue(col, out value);
        }
    }
}