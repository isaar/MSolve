using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices.Builders
{
    public class PartialMatrixRows
    {
        private readonly int numColumns;
        private readonly int numRows;
        private Dictionary<int, Dictionary<int, double>> rows;

        public PartialMatrixRows(int numRows, int numColumns)
        {
            this.numRows = numRows;
            this.numColumns = numColumns;
            this.rows = new Dictionary<int, Dictionary<int, double>>();
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
            if (rows.TryGetValue(rowIdx, out Dictionary<int, double> wholeRow)) // The row exists. Mutate it.
            {
                // If the entry did not exist, oldValue = default(double) = 0.0 and adding it should not reduce precision.
                wholeRow.TryGetValue(colIdx, out double oldValue);
                wholeRow[colIdx] = value + oldValue;
                //The Dictionary wholeRow is indexed twice in both cases. Is it possible to only index it once?
            }
            else // The row did not exist. Create a new one with the desired entry.
            {
                wholeRow = new Dictionary<int, double>();
                wholeRow.Add(colIdx, value);
                rows.Add(rowIdx, wholeRow);
            }
        }

        public SparseVector MultiplyRight(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(numColumns, vector.Length);
            var result = new SortedDictionary<int, double>();
            foreach (var wholeRow in rows)
            {
                double dot = 0.0;
                foreach (var colValPair in wholeRow.Value)
                {
                    dot += colValPair.Value * vector[colValPair.Key];
                }
                result[wholeRow.Key] = dot;
            }
            return SparseVector.CreateFromDictionary(numRows, result);
        }
    }
}
