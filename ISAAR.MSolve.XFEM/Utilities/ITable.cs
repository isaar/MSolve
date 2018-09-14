using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Utilities
{
    public interface ITable<TRow, TColumn, TValue>: IEnumerable<Tuple<TRow, TColumn, TValue>>
    {
        int EntryCount { get; }
        TValue this[TRow row, TColumn col] { get; }
        bool Contains(TRow row, TColumn col);
        IEnumerable<TColumn> GetColumnsOfRow(TRow row);
        IEnumerable<TRow> GetRows();
        IEnumerable<TValue> GetValuesOfRow(TRow row);
        bool TryGetValue(TRow row, TColumn col, out TValue value);
    }
}
