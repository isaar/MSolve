using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: Should sparsity pattern classes be together with MatrixBuilders, in a Reordering namespace in LinearAlgebra or in a 
//      Reordering namespace in Solvers?
//TODO: SuiteSparse ignores the diagonals (needs verification). Can I also ignore them when building the pattern? Update the docs
//TODO: Find a unified ordering interface and abstract the calls
namespace ISAAR.MSolve.LinearAlgebra.Reordering
{
    /// <summary>
    /// Builder for the sparsity pattern of a symmetric matrix. If the values of the matrix entries are needed and the sparsity 
    /// pattern will not change, use a DOK matrix instead. Only the super-diagonal part is stored (including the diagonal). 
    /// Efficient for outputting the pattern to column major data structures (e.g. CSC arrays).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SparsityPatternSymmetric: ISparsityPatternSymmetric
    {
        private readonly int order;
        private readonly HashSet<int>[] columns; //see performance concerns in DOKSymmetricColMajor.columns

        private SparsityPatternSymmetric(int order, HashSet<int>[] columns)
        {
            this.order = order;
            this.columns = columns;
        }

        /// <summary>
        /// See <see cref="ISparsityPattern.NumColumns"/>.
        /// </summary>
        public int NumColumns { get { return order; } }

        /// <summary>
        /// See <see cref="ISparsityPattern.NumRows"/>.
        /// </summary>
        public int NumRows { get { return order; } }

        /// <summary>
        /// The number of rows/columns of the square matrix.
        /// </summary>
        public int Order { get { return order; } }

        /// <summary>
        /// Initializes a new instance of <see cref="SparsityPatternSymmetric"/> for a 
        /// <paramref name="order"/>-by-<paramref name="order"/> symmetric matrix. Initially all entries are assumed to be 0.
        /// </summary>
        /// <param name="order">The number of rows/columns of the matrix.</param>
        public static SparsityPatternSymmetric CreateEmpty(int order)
        {
            var columns = new HashSet<int>[order];
            for (int j = 0; j < order; ++j) columns[j] = new HashSet<int>(); //Initial capacity may be optimized.
            return new SparsityPatternSymmetric(order, columns);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="SparsityPatternSymmetric"/> for an existing matrix
        /// provided as <paramref name="dense"/>.
        /// Initialy, only the entries of <paramref name="dense"/> that satisfy: Math.Abs(<paramref name="dense"/>[i, j]) &gt;
        /// <paramref name="tolerance"/> will be considered as non-zero. 
        /// </summary>
        /// <param name="dense">A matrix whose sparsity pattern will be extracted. It must be square and symmetric. Only its
        ///     upper triangle entries will be processed, thus its symmetry will not be checked explicitly.</param>
        /// <param name="tolerance">The tolerance under which an entry of <paramref name="dense"/> is considered as 0.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if <paramref name="dense"/> is not 
        ///     square.</exception>
        public static SparsityPatternSymmetric CreateFromDense(Matrix dense, double tolerance = 0.0)
        {
            Preconditions.CheckSquare(dense);
            SparsityPatternSymmetric pattern = CreateEmpty(dense.NumColumns);
            if (tolerance == 0.0)
            {
                for (int j = 0; j < dense.NumColumns; ++j)
                {
                    for (int i = 0; i <= j; ++i)
                    {
                        if (dense[i, j] != 0) pattern.columns[j].Add(i);
                    }
                }
            }
            else
            {
                for (int j = 0; j < dense.NumColumns; ++j)
                {
                    for (int i = 0; i <= j; ++i)
                    {
                        if (Math.Abs(dense[i, j]) > tolerance) pattern.columns[j].Add(i);
                    }
                }
            }
            return pattern;
        }

        /// <summary>
        /// Marks the matrix entry (<paramref name="rowIdx"/>, <paramref name="colIdx"/>) as non-zero. If it was already marked,
        /// nothing happens.
        /// </summary>
        /// <param name="rowIdx">The row index of the new non-zero entry.</param>
        /// <param name="colIdx">The column index of the new non-zero entry.</param>
        public void AddEntry(int rowIdx, int colIdx)
        {
            if (rowIdx > colIdx)
            {
                int swap = rowIdx;
                rowIdx = colIdx;
                colIdx = swap;
            }
            columns[colIdx].Add(rowIdx); //TODO: should I return whether the entry exists?
        }

        /// <summary>
        /// Like <see cref="AddEntry(int, int)"/>, but will not check or transpose the entry.
        /// </summary>
        /// <param name="rowIdx">The row index of the new non-zero entry.</param>
        /// <param name="colIdx">The column index of the new non-zero entry.</para
        public void AddEntryUpper(int rowIdx, int colIdx) => columns[colIdx].Add(rowIdx);

        /// <summary>
        /// Marks all possible matrix entries (i, j), such that i, j belong to <paramref name="indices"/>, as non-zero.
        /// </summary>
        /// <param name="indices">The row/col indices of the new non-zero entries. This array will be sorted, if it isn't 
        ///     already.</param>
        /// <param name="sorted">The new indices must be sorted before processing them. If the array <paramref name="indices"/> 
        ///     is sorted, set <paramref name="sorted"/> to true, in order to avoid resorting it. Otherwise set
        ///     <paramref name="sorted"/> to false.</param>
        public void ConnectIndices(int[] indices, bool sorted)
        {
            // Sorting the indices and processing n*(n+1)/2 entries is much faster than accessing n*n entries, transposing half
            if (!sorted) Array.Sort(indices);
            for (int j = 0; j < indices.Length; ++j)
            {
                int col = indices[j];
                for (int i = 0; i <= j; ++i) columns[col].Add(indices[i]);
            }
        }

        /// <summary>
        /// Marks all possible matrix entries (i, j), such that i, j belong to <paramref name="indices"/>, as non-zero.
        /// </summary>
        /// <param name="indices">The row/col indices of the new non-zero entries. This list will be sorted, if it isn't 
        ///     already.</param>
        /// <param name="sorted">The new indices must be sorted before processing them. If the list <paramref name="indices"/> 
        ///     is sorted, set <paramref name="sorted"/> to true, in order to avoid resorting it. Otherwise set
        ///     <paramref name="sorted"/> to false.</param>
        public void ConnectIndices(List<int> indices, bool sorted) //TODO: could I get away with OrderBy()
        {
            // Sorting the indices and processing n*(n+1)/2 entries is much faster than accessing n*n entries, transposing half
            if (!sorted) indices.Sort();
            for (int j = 0; j < indices.Count; ++j)
            {
                int col = indices[j];
                for (int i = 0; i <= j; ++i) columns[col].Add(indices[i]);
            }
        }

        /// <summary>
        /// See <see cref="ISparsityPattern.CountNonZeros"/>.
        /// </summary>
        public int CountNonZeros()
        {
            int count = 0;
            for (int j = 0; j < order; ++j)
            {
                //count += (columns[j].Count - 1) * 2 + 1; // Faster, but doesn't work if the diagonal is not present.
                foreach (int i in columns[j])
                {
                    if (i == j) ++count;
                    else count += 2;
                }
            }
            return count;
        }

        /// <summary>
        /// See <see cref="ISparsityPatternSymmetric.CountNonZeros"/>.
        /// </summary>
        public int CountNonZerosUpper()
        {
            int count = 0;
            for (int j = 0; j < order; ++j) count += columns[j].Count;
            return count;
        }

        /// <summary>
        /// See <see cref="ISparsityPattern.EnumerateNonZeros"/>.
        /// </summary>
        public IEnumerable<(int row, int col)> EnumerateNonZeros()
        {
            for (int j = 0; j < order; ++j)
            {
                foreach (var i in columns[j])
                {
                    if (i == j) yield return (i, j);
                    else //Each upper triangle entries has a corresponding lower triangle entry.
                    {
                        yield return (i, j);
                        yield return (j, i);
                    }
                }
            }
        }

        /// <summary>
        /// See <see cref="ISparsityPatternSymmetric.EnumerateNonZerosUpper"/>.
        /// </summary>
        public IEnumerable<(int row, int col)> EnumerateNonZerosUpper()
        {
            for (int j = 0; j < order; ++j)
            {
                foreach (var i in columns[j]) yield return (i, j);
            }
        }

        /// <summary>
        /// See <see cref="ISparsityPattern.IsNonZero(int, int)"/>.
        /// </summary>
        public bool IsNonZero(int rowIdx, int colIdx)
        {
            if (rowIdx > colIdx)
            {
                int swap = rowIdx;
                rowIdx = colIdx;
                colIdx = swap;
            }
            return columns[colIdx].Contains(rowIdx);
        }

        /// <summary>
        /// Assembles the indexing arrays according to the Compressed Sparse Columns format.
        /// </summary>
        /// <param name="sortRowsOfEachCol">If true, the row indices of the same column will be in ascending order. If false, 
        ///     their order will be indefinite. Sorting these indices adds an overhead, but will increase performance of the 
        ///     resulting CSC matrix's operations. It may also be required for some operations.</param>
        /// <returns>For a description of the CSC indexing arrays see <see cref="CscMatrix"/>.</returns>
        internal (int[] rowIndices, int[] colOffsets) BuildSymmetricCSCArrays(bool sortRowsOfEachCol)
        {
            int[] colOffsets = new int[order + 1];
            int nnz = 0;
            for (int j = 0; j < order; ++j)
            {
                colOffsets[j] = nnz;
                nnz += columns[j].Count;
            }
            colOffsets[order] = nnz; //The last CSC entry is nnz.

            int[] rowIndices = new int[nnz];
            double[] values = new double[nnz];
            int counter = 0;
            if (sortRowsOfEachCol)
            {
                for (int j = 0; j < order; ++j)
                {
                    //TODO: LINQ's OrderBy() might be faster and require less memory
                    foreach (var rowIdx in new SortedSet<int>(columns[j])) 
                    {
                        rowIndices[counter] = rowIdx;
                        ++counter;
                    }
                }
            }
            else
            {
                for (int j = 0; j < order; ++j)
                {
                    foreach (var rowIdx in columns[j])
                    {
                        rowIndices[counter] = rowIdx;
                        ++counter;
                    }
                }
            }

            return (rowIndices, colOffsets);
        }
    }
}
