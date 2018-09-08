using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using ISAAR.MSolve.LinearAlgebra.Exceptions;

//TODO: Should the AddSubmatrix...() methods checks that the submatrix respects the pattern? 
namespace ISAAR.MSolve.LinearAlgebra.Matrices.Builders
{
    /// <summary>
    /// Use this class for building large symmetric sparse matrices not for operations. Convert to other matrix formats once 
    /// finished and use them instead for matrix operations. The large matrices and their properties will be characterized as 
    /// "global" in this namespace. This class is meant for building symmetric global matrices in skyline storage format
    /// (see <see cref="SkylineMatrix"/>), where only entries of the upper triangle are explicitly stored.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SkylineBuilder: ISymmetricMatrixBuilder
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

        /// <summary>
        /// The number of columns of the matrix.
        /// </summary>
        public int NumColumns { get => order; }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get => order; }

        /// <summary>
        /// See <see cref="IIndexable2D.this[int, int]"/>, <see cref="IMatrixBuilder.this[int, int]"/>.
        /// </summary>
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
        /// Initializes a new instance of <see cref="SkylineBuilder"/> with the provided dimensions and sparsity bandwidth.
        /// </summary>
        /// <param name="numColumns">The number of columns in the matrix.</param>
        /// <param name="columnHeights">An array of length <paramref name="numColumns"/> containing the height of each column, 
        ///     namely the distance between the diagonal entry (exclusive) and the non-zero entry with the minimum row index 
        ///     (inclusive) of that column. The first entry must be <paramref name="columnHeights"/>[0] = 0.</param>
        /// <exception cref="ArgumentException">Thrown if the number of columns defined by <paramref name="numColumns"/>
        ///     and <paramref name="columnHeights"/> is different or if <paramref name="columnHeights"/>[0] != 0.</exception>
        public static SkylineBuilder Create(int numColumns, int[] columnHeights)
        {
            if (columnHeights.Length != numColumns) throw new ArgumentException(
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
        /// See <see cref="IMatrixBuilder.AddToEntry(int, int, double)"/>
        /// </summary>
        public void AddToEntry(int rowIdx, int colIdx, double value)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// See <see cref="ISymmetricMatrixBuilder.AddSubmatrix(IIndexable2D, IReadOnlyDictionary{int, int}, 
        /// IReadOnlyDictionary{int, int})"/>.
        /// </summary>
        /// <remarks>The caller is respondible for the global indices falling inside the skyline sparsity pattern.</remarks>
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
        /// See <see cref="ISymmetricMatrixBuilder.AddSubmatrixSymmetric(IIndexable2D, IReadOnlyDictionary{int, int})"/>.
        /// </summary>
        /// <remarks>The caller is respondible for the global indices falling inside the skyline sparsity pattern.</remarks>
        public void AddSubmatrixSymmetric(IIndexable2D subMatrix, IReadOnlyDictionary<int, int> subIndicesToGlobalIndices)
        {
            foreach (var colPair in subIndicesToGlobalIndices)
            {
                int subCol = colPair.Key;
                int globalCol = colPair.Value;
                foreach (var rowPair in subIndicesToGlobalIndices)
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
        /// See <see cref="ISymmetricMatrixBuilder.AddSubmatrixToLowerTriangle(IIndexable2D, IReadOnlyDictionary{int, int}, 
        /// IReadOnlyDictionary{int, int})"/>
        /// </summary>
        /// <remarks>The caller is respondible for the global indices falling inside the skyline sparsity pattern.</remarks>
        public void AddSubmatrixToLowerTriangle(IIndexable2D subMatrix,
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

        /// <summary>
        /// See <see cref="ISymmetricMatrixBuilder.AddSubmatrixToUpperTriangle(IIndexable2D, IReadOnlyDictionary{int, int}, 
        /// IReadOnlyDictionary{int, int})"/>
        /// </summary>
        /// <remarks>The caller is respondible for the global indices falling inside the skyline sparsity pattern.</remarks>
        public void AddSubmatrixToUpperTriangle(IIndexable2D subMatrix,
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
        /// Initializes a <see cref="SkylineMatrix"/> representation of the current matrix. This method should be 
        /// called after fully defining the matrix in <see cref="SkylineBuilder"/> format.
        /// </summary>
        /// <exception cref="EmptyMatrixBuilderException">Thrown if no non-zero entries have been defined yet.</exception>
        public SkylineMatrix BuildSkylineMatrix()
        {
            if (values.Length == 0) throw new EmptyMatrixBuilderException("Cannot build skyline arrays from a DOK with nnz = 0.");
            return SkylineMatrix.CreateFromArrays(order, values, diagOffsets, false, false);
        }

        /// <summary>
        /// See <see cref="IIndexable2D.Equals(IIndexable2D, double)"/>.
        /// </summary>
        public bool Equals(IIndexable2D other, double tolerance = 1E-13)
        {
            return BuildSkylineMatrix().Equals(other, tolerance);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void ProcessIndices(ref int rowIdx, ref int colIdx)
        { //TODO: Perhaps this should be manually inlined. Testing needed.
            if (rowIdx > colIdx)
            {
                int swap = rowIdx;
                rowIdx = colIdx;
                colIdx = swap;
            }
        }
    }
}
