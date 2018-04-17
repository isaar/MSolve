using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices.Builders
{
    /// <summary>
    /// Use this class for building a sparse matrix, e.g. <see cref="CSCMatrix"/> not for operations. Convert to other matrix 
    /// formats once finished and use them instead for matrix operations. Only the non zero entries of the upper triangle are  
    /// stored. This class is optimized for building global positive definite matrices, where there is at least 1 entry per 
    /// column.
    /// </summary>
    public class DOKColMajor: ISparseMatrix
    {
        /// <summary>
        /// See the rant in <see cref="DOKSymmetricColMajor.columns"/> about performance.
        /// </summary>
        private readonly Dictionary<int, double>[] columns;

        private DOKColMajor(int numRows, int numCols, Dictionary<int, double>[] columns)
        {
            this.columns = columns;
            this.NumRows = numRows;
            this.NumColumns = numCols;
        }

        public int NumColumns { get; }
        public int NumRows { get; }

        // Perhaps I should check indices
        public double this[int rowIdx, int colIdx]
        {
            get
            {
                if (columns[colIdx].TryGetValue(rowIdx, out double val)) return val;
                else return 0.0;
            }
            set //not thread safe
            {
                columns[colIdx][rowIdx] = value;
            }
        }

        public static DOKColMajor CreateEmpty(int numRows, int numCols)
        {
            var columns = new Dictionary<int, double>[numCols];
            for (int j = 0; j < numCols; ++j) columns[j] = new Dictionary<int, double>(); //Initial capacity may be optimized.
            return new DOKColMajor(numRows, numCols, columns);
        }

        public static DOKColMajor CreateIdentity(int order)
        {
            var columns = new Dictionary<int, double>[order];
            for (int j = 0; j < order; ++j)
            {
                var idenityCol = new Dictionary<int, double>(); //Initial capacity may be optimized.
                idenityCol[j] = 1.0;
                columns[j] = idenityCol;
            }
            return new DOKColMajor(order, order, columns);
        }

        #region global matrix building 
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
            if (columns[colIdx].TryGetValue(rowIdx, out double oldValue))
            {
                columns[colIdx][rowIdx] = value + oldValue;
            }
            else columns[colIdx][rowIdx] = value;
            //The Dictionary columns[rowIdx] is indexed twice in both cases. Is it possible to only index it once?
        }

        public void AddSubmatrix(IIndexable2D elementMatrix, 
            IReadOnlyDictionary<int, int> elementRowsToGlobalRows, IReadOnlyDictionary<int, int> elementColsToGlobalCols)
        {
            foreach (var colPair in elementColsToGlobalCols) // Col major ordering is the default
            {
                int elementCol = colPair.Key;
                int globalCol = colPair.Value;
                foreach (var rowPair in elementRowsToGlobalRows)
                {
                    int elementRow = rowPair.Key;
                    int globalRow = rowPair.Value;
                    double elementValue = elementMatrix[elementRow, elementCol];
                    if (columns[globalCol].TryGetValue(globalRow, out double oldGlobalValue))
                    {
                        columns[globalCol][globalRow] = elementValue + oldGlobalValue;
                    }
                    else columns[globalCol][globalRow] = elementValue;
                }
            }
        }

        /// <summary>
        /// Use this method if 1) both the global DOK and the element matrix are symmetric and 2) the rows and columns correspond
        /// to the same degrees of freedom. The caller is responsible for making sure that both matrices are symmetric and that 
        /// the dimensions and dofs of the element matrix and dof mappings match.
        /// </summary>
        /// <param name="elementMatrix">The element submatrix, entries of which will be added to the global DOK. It must be 
        ///     symmetric and its <see cref="IIndexable2D.NumColumns"/> = <see cref="IIndexable2D.NumRows"/> must be equal to
        ///     elemenDofs.Length = globalDofs.Length.</param>
        /// <param name="elementDofs">The entries in the element matrix to be added to the global matrix. Specificaly, pairs of 
        ///     (elementDofs[i], elementDofs[j]) will be added to (globalDofs[i], globalDofs[j]).</param>
        /// <param name="globalDofs">The entries in the global matrix where element matrix entries will be added to. Specificaly,
        ///     pairs of (elementDofs[i], elementDofs[j]) will be added to (globalDofs[i], globalDofs[j]).</param>
        public void AddSubmatrixSymmetric(IIndexable2D elementMatrix, int[] elementDofs, int[] globalDofs)
        {
            int n = elementDofs.Length;
            for (int j = 0; j < n; ++j)
            {
                int elementCol = elementDofs[j];
                int globalCol = globalDofs[j];

                //Diagonal entry
                if (columns[globalCol].TryGetValue(globalCol, out double oldGlobalDiagValue))
                {
                    columns[globalCol][globalCol] = elementMatrix[elementCol, elementCol] + oldGlobalDiagValue;
                }
                else columns[globalCol][globalCol] = elementMatrix[elementCol, elementCol];

                //Non diagonal entries
                for (int i = 0; i < j; ++i)
                {
                    int elementRow = elementDofs[j];
                    int globalRow = globalDofs[j];
                    double newGlobalValue = elementMatrix[elementRow, elementCol];
                    // Only check the upper triangle. If the DOK matrix is not symmetric, this will cause errors
                    if (columns[globalCol].TryGetValue(globalRow, out double oldGlobalValue)) 
                    {
                        newGlobalValue += oldGlobalValue;
                    }
                    columns[globalCol][globalRow] = newGlobalValue;
                    columns[globalRow][globalCol] = newGlobalValue;
                }
            }
        }
        #endregion

        /// <summary>
        /// Creates the CSC arrays.
        /// </summary>
        /// <param name="sortRowsOfEachCol">True to sort the row indices of the CSC matrix between colOffsets[j] and 
        ///     colOffsets[j+1] in ascending order. False to leave them unordered. Ordered rows might result in better  
        ///     performance during multiplications or they might be required by 3rd party libraries. Leaving them unordered will 
        ///     be faster during creation of the CSC matrix.</param>
        /// <returns></returns>
        public (double[] values, int[] rowIndices, int[] columnOffsets) BuildCSCArrays(bool sortRowsOfEachCol)
        {
            int[] colOffsets = new int[NumColumns + 1];
            int nnz = 0;
            for (int j = 0; j < NumColumns; ++j)
            {
                colOffsets[j] = nnz;
                nnz += columns[j].Count;
            }
            colOffsets[NumColumns] = nnz; //The last CSC entry is nnz.

            int[] rowIndices = new int[nnz];
            double[] values = new double[nnz];
            int counter = 0;
            for (int j = 0; j < NumColumns; ++j)
            {
                if (sortRowsOfEachCol)
                {
                    foreach (var rowVal in columns[j].OrderBy(pair => pair.Key))
                    {
                        rowIndices[counter] = rowVal.Key;
                        values[counter] = rowVal.Value;
                        ++counter;
                    }
                }
                else
                {
                    foreach (var rowVal in columns[j])
                    {
                        rowIndices[counter] = rowVal.Key;
                        values[counter] = rowVal.Value;
                        ++counter;
                    }
                }
            }

            return (values, rowIndices, colOffsets);
        }

        /// <summary>
        /// Creates a CSC matrix.
        /// </summary>
        /// <param name="sortRowsOfEachCol">True to sort the row indices of the CSC matrix between colOffsets[j] and 
        ///     colOffsets[j+1] in ascending order. False to leave them unordered. Ordered rows might result in better  
        ///     performance during multiplications or they might be required by 3rd party libraries. Leaving them unordered will 
        ///     be faster during creation of the CSC matrix.</param>
        /// <returns></returns>
        public CSCMatrix BuildCSCMatrix(bool sortRowsOfEachCol)
        {
            (double[] values, int[] rowIndices, int[] colOffsets) = BuildCSCArrays(sortRowsOfEachCol);
            return CSCMatrix.CreateFromArrays(NumRows, NumColumns, values, rowIndices, colOffsets, false);
        }

        public int CountNonZeros()
        {
            int count = 0;
            for (int j = 0; j < NumColumns; ++j) count += columns[j].Count;
            return count;
        }

        public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
        {
            for (int j = 0; j < NumColumns; ++j)
            {
                foreach (var rowVal in columns[j])
                {
                    yield return (rowVal.Key, j, rowVal.Value);
                }
            }
        }

        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
        {
            return DenseStrategies.AreEqual(this, other, tolerance);
        }

        public SparseFormat GetSparseFormat()
        {
            (double[] values, int[] rowIndices, int[] colOffsets) = BuildCSCArrays(false);
            var format = new SparseFormat();
            format.RawValuesTitle = "Values";
            format.RawValuesArray = values;
            format.RawIndexArrays.Add("Row indices", rowIndices);
            format.RawIndexArrays.Add("Column offsets", colOffsets);
            return format;
        }
    }
}
