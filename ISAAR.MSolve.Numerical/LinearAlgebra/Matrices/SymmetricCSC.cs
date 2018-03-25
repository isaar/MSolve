using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    /// <summary>
    /// Compressed Sparse Columns format.
    /// </summary>
    public class SymmetricCSC: IIndexable2D
    {
        /// <summary>
        /// values.Length = number of non zeros in upper triangle.
        /// </summary>
        private readonly double[] values;

        /// <summary>
        /// rowIndices.Length = number of non zeros in upper triangle.
        /// </summary>
        private readonly int[] rowIndices;

        /// <summary>
        /// colOffsets.Length = number of rows/columns +1. 
        /// colOffsets[colOffsets.Length-1] = number of non zeros in upper triangle.
        /// </summary>
        private readonly int[] colOffsets;

        /// <summary>
        /// Does not copy anything
        /// </summary>
        /// <param name="values"></param>
        /// <param name="rowIndices"></param>
        /// <param name="colOffsets"></param>
        /// <param name="checkInput">Set to true to perform basic dimension checks.</param>
        public SymmetricCSC(double[] values, int[] rowIndices, int[] colOffsets, bool checkInput)
        {
            int nnz = colOffsets[colOffsets.Length - 1];
            if (checkInput)
            {
                if ((nnz != values.Length) || (nnz != rowIndices.Length))
                {
                    throw new ArgumentException("Mismatch in dimensions of the CSC arrays. Check that colOffsets.Length ="
                        + " number of rows/columns + 1 and that colOffsets[colOffsets.Length-1] = values.Length ="
                        + " rowIndices.Length = number of non zeros in upper triangle");
                }
            }
            this.values = values;
            this.rowIndices = rowIndices;
            this.colOffsets = colOffsets;
            this.NumColumns = colOffsets.Length - 1;
            this.NumRows = colOffsets.Length - 1;
            this.NumNonZerosUpper = nnz;
        }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// Only structural non zeros
        /// </summary>
        /// <returns></returns>
        public int NumNonZeros { get { throw new NotImplementedException(); } }

        public int NumNonZerosUpper { get; }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get; }

        public double this[int rowIdx, int colIdx]
        {
            get
            {
                ProcessIndices(ref rowIdx, ref colIdx);
                int colStart = colOffsets[colIdx];
                int colEnd = colOffsets[colIdx+1];
                for (int i = colStart; i < colEnd; ++i) //Only scan the nnz entries for the given column
                {
                    if (rowIndices[i] == rowIdx) return values[i];
                }
                return 0.0; //The requensted entry is not stored.
            }
        }

        public SuiteSparseCholesky FactorCholesky()
        {
            return SuiteSparseCholesky.CalcFactorization(NumColumns, NumNonZerosUpper, values, rowIndices, colOffsets);
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
