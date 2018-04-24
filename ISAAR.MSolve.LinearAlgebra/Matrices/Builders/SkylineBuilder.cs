using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Exceptions;

//TODO: Should the AddSubmatrix...() methods checks that the submatrix respects the pattern? 
namespace ISAAR.MSolve.LinearAlgebra.Matrices.Builders
{
    /// <summary>
    /// For symmetric and column major skyline matrices.
    /// </summary>
    public class SkylineBuilder: IIndexable2D
    {
        private readonly int order;
        private readonly int numNonZeros;
        private readonly double[] values;
        private readonly int[] diagOffsets;

        private SkylineBuilder(int order, int numNonZeros, double[] values, int[] diagOffsets)
        {
            this.order = order;
            this.numNonZeros = numNonZeros;
            this.values = values;
            this.diagOffsets = diagOffsets;
        }

        public int NumColumns { get { return order; } }
        public int NumRows { get { return order; } }

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

        /// <summary>
        /// 
        /// </summary>
        /// <param name="numColumns">The number of columns in the matrix.</param>
        /// <param name="columnHeights">An array of length <paramref name="numColumns"/> containing the height of each column, 
        ///     namely the distance between the diagonal entry (exclusive) and the non-zero entry with the minimum row index 
        ///     (inclusive) of that column. WARNING: As said before, these heights do not count the diagonal entry.</param>
        /// <returns></returns>
        public static SkylineBuilder Create(int numColumns, int[] columnHeights)
        {
            if (columnHeights.Length != numColumns) throw new NonMatchingDimensionsException(
                $"There are {columnHeights} column heights, but {numColumns} columns.");
            if (columnHeights[0] != 0) throw new ArgumentException(
                $"The first column's height must be 0, but was {columnHeights[0]}.");

            // Indexer
            int[] diagOffsets = new int[columnHeights.Length + 1]; // 1 extra ghost entry to facilitate some operations.
            diagOffsets[0] = 0; // The diagonal entries of the first 2 rows always have the same positions
            diagOffsets[1] = 1;
            for (int i = 1; i < columnHeights.Length; ++i)
            {
                diagOffsets[i + 1] = diagOffsets[i] + columnHeights[i] + 1; //+1 because the heights don't include the diagonal
            }

            // Initialize values to 0
            int nnz = diagOffsets[diagOffsets.Length - 1];
            double[] values = new double[nnz];

            return new SkylineBuilder(numColumns,nnz, values, diagOffsets);
        }

        /// <summary>
        /// Maps <paramref name="subMatrix"/> to "global" rows and columns of the underlying skyline matrix and then adds 
        /// the entries of <paramref name="subMatrix"/> to these global indices. The caller is respondible for the global 
        /// indices falling inside the skyline sparsity pattern. WARNING: this method adds the whole <paramref name="subMatrix"/>.
        /// If you only want to add the upper triangular part of it, use 
        /// <see cref="AddSubmatrixSymmetric(IIndexable2D, IReadOnlyDictionary{int, int})"/> instead.
        /// </summary>
        /// <param name="subMatrix"></param>
        /// <param name="subRowsToGlobalRows"></param>
        /// <param name="subColsToGlobalCols"></param>
        public void AddSubmatrix(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols)
        {
            foreach (var colPair in subColsToGlobalCols)
            {
                int subCol = colPair.Key;
                foreach (var rowPair in subRowsToGlobalRows)
                {
                    int globalRow, globalCol; // Postcondition: globalCol >= globalRow, thus height >=0
                    if (rowPair.Value <= colPair.Value)
                    {
                        globalRow = rowPair.Value;
                        globalCol = colPair.Value;
                    }
                    else
                    {
                        globalRow = colPair.Value;
                        globalCol = rowPair.Value;
                    }

                    int offset = diagOffsets[globalCol] + globalCol - globalRow;
                    values[offset] += subMatrix[rowPair.Key, subCol];
                }
            }
        }

        /// <summary>
        /// Maps <paramref name="subMatrix"/> to "global" rows and columns of the underlying skyline matrix and then adds 
        /// the entries of <paramref name="subMatrix"/> to these global indices. The caller is respondible for the global 
        /// indices falling inside the skyline sparsity pattern. Optimized version of 
        /// <see cref="AddSubmatrix(IIndexable2D, IReadOnlyDictionary{int, int}, IReadOnlyDictionary{int, int})"/>, in case the
        /// caller is sure that all global indices will be super diagonal.
        /// </summary>
        public void AddSubmatrixAboveDiagonal(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols)
        {
            foreach (var colPair in subColsToGlobalCols) // Col major ordering
            {
                int subCol = colPair.Key;
                int globalCol = colPair.Value;
                foreach (var rowPair in subRowsToGlobalRows)
                {
                    int globalRow = rowPair.Value;

                    // The assumption is: globalCol >= globalRow, thus height >=0
                    int offset = diagOffsets[globalCol] + globalCol - globalRow;
                    values[offset] += subMatrix[rowPair.Key, subCol];
                }
            }
        }

        /// <summary>
        /// Maps <paramref name="subMatrix"/> to "global" rows and columns of the underlying skyline matrix and then adds the
        /// entries of <paramref name="subMatrix"/> to these global indices. The caller is respondible for the global indices
        /// falling inside the skyline sparsity pattern. Use this method if you only want to add the upper triangular part of the
        /// submatrix. Otherwise <see cref="AddSubmatrix(IIndexable2D, IReadOnlyDictionary{int, int}, IReadOnlyDictionary{int, int})"/>
        /// will add the shole submatrix, resulting in the off-diagonal entries being added twice.
        /// </summary>
        public void AddSubmatrixSymmetric(IIndexable2D subMatrix, IReadOnlyDictionary<int, int> subDOFsToGlobalDOFs)
        {
            foreach (var colPair in subDOFsToGlobalDOFs)
            {
                int subCol = colPair.Key;
                int globalCol = colPair.Value;
                foreach (var rowPair in subDOFsToGlobalDOFs)
                {
                    int globalRow = rowPair.Value;

                    int height = globalCol - globalRow;
                    if (height >= 0) // Only process the superdiagonal entries
                    {
                        int offset = diagOffsets[globalCol] + globalCol - globalRow;
                        values[offset] += subMatrix[rowPair.Key, subCol];
                    }
                }
            }
        }

        /// <summary>
        /// Maps <paramref name="subMatrix"/> to "global" rows and columns of the underlying skyline matrix and then adds 
        /// the entries of <paramref name="subMatrix"/> to these global indices. The caller is respondible for the global 
        /// indices falling inside the skyline sparsity pattern. Optimized version of 
        /// <see cref="AddSubmatrix(IIndexable2D, IReadOnlyDictionary{int, int}, IReadOnlyDictionary{int, int})"/>, in case the
        /// caller is sure that all global indices will be super diagonal.
        /// </summary>
        public void AddSubmatrixUnderDiagonal(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols)
        {
            foreach (var rowPair in subRowsToGlobalRows) // Transpose(Col major ordering) = Row major ordering
            {
                int subRow = rowPair.Key;
                int globalRow = rowPair.Value;
                foreach (var colPair in subColsToGlobalCols)
                {
                    int globalCol = colPair.Value;

                    // The assumption is: globalCol <= globalRow, thus height <=0
                    // Transpose the access on the global matrix, instead of swapping the indices. 
                    int offset = diagOffsets[globalRow] + globalRow - globalCol;
                    values[offset] += subMatrix[subRow, colPair.Key];
                }
            }
        }

        public SkylineMatrix BuildSkylineMatrix()
        {
            return SkylineMatrix.CreateFromArrays(order, values, diagOffsets, false, false);
        }

        public bool Equals(IIndexable2D other, double tolerance = 1E-13)
        {
            return BuildSkylineMatrix().Equals(other, tolerance);
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
