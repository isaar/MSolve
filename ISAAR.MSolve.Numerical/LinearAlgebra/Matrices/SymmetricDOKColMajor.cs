using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    /// <summary>
    /// Use this class for building a symmetric sparse matrix, not for operations. Convert to other matrix formats once finished 
    /// and use them instead for matrix operations. Only the non zero entries of the upper triangle are stored.
    /// </summary>
    public class SymmetricDOKColMajor : IIndexable2D
    {
        /// <summary>
        /// An array of dictionaries is more efficent and perhaps easier to work with than a dictionary of dictionaries. There 
        /// is usually at least 1 non zero entry in each column. Otherwise this data structure wastes space, but if there were  
        /// many empty rows, perhaps another matrix format is more appropriate.
        /// To get the value: data[colIdx][rowIdx] = value. 
        /// To get the row-value subdictionary: data[colIdx] = SortedDictionary[int, double]
        /// TODO: Perhaps a Dictionary should be used instead of SortedDictionary and only sort each column independently before 
        /// converting. Dictionary is much more efficient for insertion, which will usually happen more than once for each entry. 
        /// </summary>
        private readonly Dictionary<int, double>[] data;

        public SymmetricDOKColMajor(int numCols)
        {
            this.data = new Dictionary<int, double>[numCols];
            for (int j = 0; j < numCols; ++j) this.data[j] = new Dictionary<int, double>(); // Initial capacity may be optimized.
            this.NumColumns = numCols;
        }

        public int NumColumns { get; }
        public int NumRows { get { return NumColumns; } }

        public double this[int rowIdx, int colIdx]
        {
            get
            {
                ProcessIndices(ref rowIdx, ref colIdx);
                if (data[colIdx].TryGetValue(rowIdx, out double val)) return val;
                else return 0.0;
            }
            set //not thread safe
            {
                ProcessIndices(ref rowIdx, ref colIdx);
                data[colIdx][rowIdx] = value;
            }
        }

        /// <summary>
        /// If the entry already exists: the new value is added to the existing one: this[rowIdx, colIdx] += value. 
        /// Otherwise the entry is inserted with the new value: this[rowIdx, colIdx] = value.
        /// Use this method instead of this[rowIdx, colIdx] += value, as it is an optimized version.
        /// </summary>
        /// <param name="rowIdx"></param>
        /// <param name="colIdx"></param>
        /// <param name="value"></param>
        public void AddToEntry(int rowIdx, int colIdx, double value)
        {
            ProcessIndices(ref rowIdx, ref colIdx);
            if (data[colIdx].TryGetValue(rowIdx, out double oldValue))
            {
                data[colIdx][rowIdx] = value + oldValue;
            }
            else data[colIdx][rowIdx] = value;
            //The Dictionary data[rowIdx] is indexed twice in both cases. Is it possible to only index it once?
        }

        public SymmetricCSC ToSymmetricCSC()
        {
            int[] colOffsets = new int[NumColumns + 1];
            int nnz = 0;
            for (int j = 0; j < NumColumns; ++j)
            {
                colOffsets[j] = nnz;
                nnz += data[j].Count;
            }
            colOffsets[NumColumns] = nnz; //The last CSC entry is nnz.

            int[] rowIndices = new int[nnz];
            double[] values = new double[nnz];
            int counter = 0;
            for (int j = 0; j < NumColumns; ++j)
            {
                foreach (var rowVal in data[j])
                {
                    rowIndices[counter] = rowVal.Key;
                    values[counter] = rowVal.Value;
                    ++counter;
                }
            }

            return new SymmetricCSC(values, rowIndices, colOffsets, false);
        }

        /// <summary>
        /// Perhaps this should be manually inlined. Testing needed.
        /// </summary>
        /// <param name="rowIdx"></param>
        /// <param name="colIdx"></param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void ProcessIndices(ref int rowIdx, ref int colIdx)
        {
            if (rowIdx > colIdx)
            {
                int swap = rowIdx;
                rowIdx = colIdx;
                colIdx = swap;
            }
        }
    }
}
