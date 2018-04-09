using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    public class SkylineMatrix: IIndexable2D, ISparseMatrix
    {
        /// <summary>
        /// Contains the non zero superdiagonal entries of the matrix in column major order, starting from the diagonal and going
        /// upwards. Its length is nnz.
        /// </summary>
        private readonly double[] values;

        /// <summary>
        /// Contains the indices into values of the diagonal entries of the matrix. Its length = order + 1, with the last entry
        /// being equal to nnz.
        /// </summary>
        private readonly int[] diagOffsets;

        private SkylineMatrix(int order, double[] values, int[] diagOffsets)
        {
            this.values = values;
            this.diagOffsets = diagOffsets;
            this.NumColumns = order;
        }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get { return NumColumns; } }

        // TODO: should I add index bound checking?
        public double this[int rowIdx, int colIdx]
        {
            get
            {
                ProcessIndices(ref rowIdx, ref colIdx);
                int entryHeight = colIdx - rowIdx; // excluding diagonal
                int maxColumnHeight = diagOffsets[colIdx + 1] - diagOffsets[colIdx] - 1; // excluding diagonal
                if (entryHeight > maxColumnHeight) return 0.0; // outside stored non zero pattern
                else return values[diagOffsets[colIdx] + entryHeight];
            }
            set
            {
                ProcessIndices(ref rowIdx, ref colIdx);
                int entryHeight = colIdx - rowIdx; // excluding diagonal
                int maxColumnHeight = diagOffsets[colIdx + 1] - diagOffsets[colIdx] - 1; // excluding diagonal
                if (entryHeight > maxColumnHeight)
                {
                    throw new SparsityPatternModifiedException($"In column {colIdx} only rows [{maxColumnHeight}, {colIdx}]"
                        + $" can be changed, but you are trying to set entry ({rowIdx}, {colIdx})");
                }
                else values[diagOffsets[colIdx] + entryHeight] = value;
            }
        }

        public static SkylineMatrix CreateFromArrays(int order, double[] values, int[] diagOffsets, bool checkInput)
        {
            if (checkInput)
            {
                if (diagOffsets.Length != order + 1)
                {
                    throw new ArgumentException("The length of the Skyline diagonal offsets array must be equal to the number of"
                        + " rows/columns + 1, but was " + diagOffsets.Length);
                }
                if (diagOffsets[diagOffsets.Length - 1] != values.Length)
                {
                    throw new ArgumentException("The last entry of the Skyline diagonal offsets array must be equal to the number"
                        + " of non zero entries, but was " + diagOffsets[diagOffsets.Length - 1]);
                }
            }
            return new SkylineMatrix(order, values, diagOffsets);
        }

        public static SkylineMatrix CreateZeroWithPattern(int order, int[] diagOffsets, bool checkInput)
        {
            if (checkInput)
            {
                if (diagOffsets.Length != order + 1)
                {
                    throw new ArgumentException("The length of the Skyline diagonal offsets array must be equal to the number of"
                        + " rows/columns + 1, but was " + diagOffsets.Length);
                }
            }
            int nnz = diagOffsets[diagOffsets.Length] - 1;
            return new SkylineMatrix(order, new double[nnz], diagOffsets);
        }

        public int CountNonZeros()
        {
            return values.Length;
        }

        public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
        {
            for (int j = 0; j < NumColumns; ++j)
            {
                int colOffset = diagOffsets[j];
                int columnTop = j - diagOffsets[j + 1] + colOffset + 1;
                yield return (j, j, values[colOffset]); // diagonal entry
                for (int i = columnTop; i < j; ++i)
                {
                    double value = colOffset + j - i;
                    yield return (i, j, value);
                    yield return (j, i, value);
                }
            }
        }

        public bool Equals(IIndexable2D other, double tolerance = 1E-13)
        {
            if ((this.NumRows != other.NumRows) || (this.NumColumns != other.NumColumns)) return false;
            var comparer = new ValueComparer(1e-13);
            for (int j = 0; j < NumColumns; ++j)
            {
                int colOffset = diagOffsets[j];
                int columnTop = j - diagOffsets[j+1] + colOffset + 1;
                for (int i = 0; i < columnTop; ++i) // zero entries above stored column
                {
                    if (!( comparer.AreEqual(0.0, other[i, j]) && comparer.AreEqual(0.0, other[j, i]) )) return false;
                }
                for (int i = columnTop; i < j; ++i) // non zero entries of column, excluding diafonal
                {
                    double value = colOffset + j - i;
                    if (!(comparer.AreEqual(value, other[i, j]) && comparer.AreEqual(value, other[j, i]))) return false;
                }
                if (!comparer.AreEqual(values[colOffset], other[j, j])) return false; // non zero diagonal entry
            }
            return true; // At this point all entries have been checked and are equal
        }

        public SparseFormat GetSparseFormat()
        {
            var format = new SparseFormat();
            format.RawValuesTitle = "Values";
            format.RawValuesArray = values;
            format.RawIndexArrays.Add("Diagonal offsets", diagOffsets);
            return format;
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
