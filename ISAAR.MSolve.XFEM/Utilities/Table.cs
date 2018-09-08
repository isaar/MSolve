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
    // TODO: If we know that the rows are ALL the nodes of the model and they are numberd contiguously, we can implement the rows 
    // as an array instead of Dictionary, in order to speed up things.
    public class Table<TRow, TColumn, TValue> : ITable<TRow, TColumn, TValue>
    {
        protected readonly Dictionary<TRow, Dictionary<TColumn, TValue>> data;

        public Table()
        {
            this.data = new Dictionary<TRow, Dictionary<TColumn, TValue>>();
        }

        protected Table(Dictionary<TRow, Dictionary<TColumn, TValue>> data)
        {
            this.data = data;
        }

        public int EntryCount //TODO: perhaps this should be cached somehow
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
                bool rowExists = data.TryGetValue(row, out Dictionary<TColumn, TValue> wholeRow);
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
            bool containsRow = data.TryGetValue(row, out Dictionary<TColumn, TValue> wholeRow);
            if (!containsRow) return false; 
            else return wholeRow.ContainsKey(col);
        }

        //TODO: use named tuple
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
            bool containsRow = data.TryGetValue(row, out Dictionary<TColumn, TValue> wholeRow);
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