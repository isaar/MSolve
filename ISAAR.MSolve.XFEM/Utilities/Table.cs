using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Utilities
{
    // TODO: Add slicing features
    // TODO: Perhaps a Table.Builder would be better, since EntryCount can be O(1)
    class Table<TRow, TColumn, TValue> : ITable<TRow, TColumn, TValue>
    {
        protected readonly Dictionary<TRow, Dictionary<TColumn, TValue>> data;

        public Table()
        {
            this.data = new Dictionary<TRow, Dictionary<TColumn, TValue>>();
        }

        public int EntryCount
        {
            get
            {
                int count = 0;
                foreach (var wholeRow in data)
                {
                    count += wholeRow.Value.Count;
                }
                return count;
            }
        }

        public TValue this[TRow row, TColumn col]
        {
            get { return data[row][col]; }
            set
            {
                Dictionary<TColumn, TValue> wholeRow;
                bool rowExists = data.TryGetValue(row, out wholeRow);
                if (!rowExists)
                {
                    wholeRow = new Dictionary<TColumn, TValue>();
                    data.Add(row, wholeRow);
                }
                if (wholeRow.ContainsKey(col))
                {
                    throw new ArgumentException("The entry (row, column) = (" + row + ", "
                        + col + ") already exists.");
                }
                else wholeRow.Add(col, value);
            }
        }

        public bool Contains(TRow row, TColumn col)
        {
            Dictionary<TColumn, TValue> wholeRow;
            bool containsRow = data.TryGetValue(row, out wholeRow);
            if (!containsRow) return false; 
            else return wholeRow.ContainsKey(col);
        }

        public IEnumerator<Tuple<TRow, TColumn, TValue>> GetEnumerator()
        {
            foreach (var wholeRow in data)
            {
                foreach (var colValPair in wholeRow.Value)
                {
                    yield return new Tuple<TRow, TColumn, TValue>(wholeRow.Key, colValPair.Key, colValPair.Value);
                }
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
        
        public IEnumerable<TColumn> GetColumnsOfRow(TRow row)
        {
            return data[row].Keys;
        }

        public IEnumerable<TRow> GetRows()
        {
            return data.Keys;
        }

        public IEnumerable<TValue> GetValuesOfRow(TRow row)
        {
            return data[row].Values;
        }

        public bool TryGetValue(TRow row, TColumn col, out TValue value)
        {
            Dictionary<TColumn, TValue> wholeRow;
            bool containsRow = data.TryGetValue(row, out wholeRow);
            if (!containsRow)
            {
                value = default(TValue);
                return false;
            }
            else return wholeRow.TryGetValue(col, out value);
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
    }
}