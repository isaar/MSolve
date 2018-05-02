using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: Should sparsity pattern classes be together with MatrixBuilders, in a Reordering namespace in LinearAlgebra or in a 
//      Reordering namespace in Solvers?
//TODO: SuiteSparse ignores the diagonals (needs verification). Can I also ignore them when building the pattern? Update the docs
namespace ISAAR.MSolve.LinearAlgebra.Reordering
{
    /// <summary>
    /// Sparsity pattern for symmetric matrices. Provides methods for building the sparsity pattern efficiently and easily. If 
    /// the whole matrix is needed and the sparsity pattern will not change, use a DOK matrix instead. 
    /// Only the super-diagonal part is stored (including the diagonal). 
    /// Efficient for outputting the pattern to column major data structures (e.g. CSC arrays).
    /// </summary>
    public class SparsityPatternSymmetricColMajor: ISparsityPatternSymmetric
    {
        private readonly int order;
        private readonly HashSet<int>[] columns; //see performance concerns in DOKSymmetricColMajor.columns

        private SparsityPatternSymmetricColMajor(int order, HashSet<int>[] columns)
        {
            this.order = order;
            this.columns = columns;
        }

        public int NumColumns { get { return order; } }
        public int NumRows { get { return order; } }
        public int Order { get { return order; } }

        public static SparsityPatternSymmetricColMajor CreateEmpty(int order)
        {
            var columns = new HashSet<int>[order];
            for (int j = 0; j < order; ++j) columns[j] = new HashSet<int>(); //Initial capacity may be optimized.
            return new SparsityPatternSymmetricColMajor(order, columns);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="dense">Must be square and symmetric. The upper triangle will be processed without checking the 
        ///     symmetricity.</param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        public static SparsityPatternSymmetricColMajor CreateFromDense(Matrix dense, double tolerance = 0.0)
        {
            Preconditions.CheckSquare(dense);
            SparsityPatternSymmetricColMajor pattern = CreateEmpty(dense.NumColumns);
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
        /// Adds a new matrix entry (<paramref name="rowIdx"/>, <paramref name="colIdx"/>). If it already existed, nothing happens.
        /// </summary>
        /// <param name="rowIdx">The row index of the new entry.</param>
        /// <param name="colIdx">The column index of the new entry.</param>
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
        /// Adds all possible matrix entries (i, j), such that i, j belong to <paramref name="indices"/>.
        /// </summary>
        /// <param name="indices">The row/col indices of the entries to be added. The array will be sorted, if it isn't 
        ///     already.</param>
        /// <param name="sorted">True if <paramref name="indices"/> are sorted. False otherwise.</param>
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

        public int CountNonZerosUpper()
        {
            int count = 0;
            for (int j = 0; j < order; ++j) count += columns[j].Count;
            return count;
        }

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

        public IEnumerable<(int row, int col)> EnumerateNonZerosUpper()
        {
            for (int j = 0; j < order; ++j)
            {
                foreach (var i in columns[j]) yield return (i, j);
            }
        }

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
        /// Returns a permutation vector using the provided ordering algorithm.
        /// </summary>
        /// <param name="orderingAlgorithm"></param>
        /// <returns></returns>
        public int[] Reorder(OrderingAMD orderingAlgorithm)
        {
            (int[] rowIndices, int[] colOffsets) = BuildSymmetricCSCArrays(true);
            return orderingAlgorithm.FindPermutation(order, rowIndices.Length, rowIndices, colOffsets);
        }

        private (int[] rowIndices, int[] colOffsets) BuildSymmetricCSCArrays(bool sortRowsOfEachCol)
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
