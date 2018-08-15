using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices.Builders
{
    /// <summary>
    /// Facilitates the construction of a square matrix when only some of its columns are needed explicitly.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PartialMatrixColumns
    {
        private readonly Dictionary<int, Dictionary<int, double>> columns;
        
        /// <summary>
        /// Initializes a new <see cref="PartialMatrixColumns"/> instance to reperesent the relevant columns of a matrix 
        /// with the specified dimensions.
        /// </summary>
        /// <param name="order">The number of rows/columns of the whole square matrix.</param>
        public PartialMatrixColumns(int order)
        {
            this.NumRows = order; ;
            this.NumColumns = order;
            this.columns = new Dictionary<int, Dictionary<int, double>>();
        }

        /// <summary>
        /// The number of columns of the whole matrix.
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of rows of the whole matrix.
        /// </summary>
        public int NumRows { get; }

        /// <summary>
        /// Adds the provided <paramref name="value"/> to the entry (<paramref name="rowIdx"/>, <paramref name="colIdx"/>). 
        /// </summary>
        /// <param name="rowIdx">The row index of the entry to modify. Constraints: 
        ///     0 &lt;= <paramref name="rowIdx"/> &lt; this.<see cref="IIndexable2D.NumRows"/>.</param>
        /// <param name="colIdx">The column index of the entry to modify. Constraints: 
        ///     0 &lt;= <paramref name="rowIdx"/> &lt; this.<see cref="IIndexable2D.NumColumns"/>.</param>
        /// <param name="value">The value that will be added to the entry (<paramref name="colIdx"/>, <paramref name="colIdx"/>).
        ///     </param>
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
        /// Returns the column with index = <paramref name="colIdx"/> as a vector. However, the entries with row index that 
        /// belongs in <paramref name="tabooRows"/> will be set to 0. More accurately, they will not be included in the 
        /// sparsity pattern of the returned <see cref="SparseVector"/>. Note that the length of the returned vector is 
        /// equal to this.<see cref="NumRows"/>.
        /// </summary>
        /// <param name="colIdx">The index of the column to return. Constraints: Column <paramref name="colIdx"/> must be stored
        ///     and 0 &lt;= <paramref name="colIdx"/> &lt; this.<see cref="NumColumns"/>.</param>
        /// <param name="tabooRows">The entries of the returned column vector at the indices specified by 
        ///     <paramref name="tabooRows"/> will be equal to 0. Constraints: foreach rowIdx in <paramref name="tabooRows"/>:
        ///     0 &lt;= rowIdx &lt; this.<see cref="NumRows"/>.</param>
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="colIdx"/> or <paramref name="tabooRows"/> 
        ///     violate the described constraints.</exception>
        public SparseVector GetColumnWithoutRows(int colIdx, HashSet<int> tabooRows)
        {
            bool exists = columns.TryGetValue(colIdx, out Dictionary<int, double> wholeColumn); //TODO: Should I just let the Dictionary indexer throw?
            if (!exists) throw new IndexOutOfRangeException($"Column {colIdx} is not stored.");
            var result = new SortedDictionary<int, double>();
            foreach (var rowVal in columns[colIdx])
            {
                if (!tabooRows.Contains(rowVal.Key)) result.Add(rowVal.Key, rowVal.Value);
            }
            return SparseVector.CreateFromDictionary(NumRows, result);
        }

        /// <summary>
        /// Returns the column with index = <paramref name="colIdx"/> as a vector. However, only the entries with row index that 
        /// belongs in <paramref name="wantedRows"/> will be copied. The rest will be 0 and not stored explicitly. Note that 
        /// the length of the returned vector is equal to this.<see cref="NumRows"/>.
        /// </summary>
        /// <param name="colIdx">The index of the column to return. Constraints: Column <paramref name="colIdx"/> must be stored
        ///     and 0 &lt;= <paramref name="colIdx"/> &lt; this.<see cref="NumColumns"/>.</param>
        /// <param name="wantedRows">The entries of the column <paramref name="colIdx"/> that will be copied to the returned 
        ///     vector. Constraints: foreach rowIdx in <paramref name="wantedRows"/>:
        ///     0 &lt;= rowIdx &lt; this.<see cref="NumRows"/>.</param>
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="colIdx"/> or <paramref name="wantedRows"/> 
        ///     violate the described constraints.</exception>
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
            return SparseVector.CreateFromDictionary(NumRows, result);
        }

        /// <summary>
        /// Sets the column with index <paramref name="colIdx"/> to have 1.0 at the main diagonal entry and 0.0 everywhere else.
        /// </summary>
        /// <param name="colIdx">The index of the column to modify. Constraints:
        ///     and 0 &lt;= <paramref name="colIdx"/> &lt; this.<see cref="NumColumns"/>.</param>
        public void SetColumnToIdentity(int colIdx)
        { //TODO: this is probably not needed as the caller could approach identity columns more efficiently
            var identityCol = new Dictionary<int, double>();
            identityCol[colIdx] = 1.0;
            columns[colIdx] = identityCol;
        }
    }
}
