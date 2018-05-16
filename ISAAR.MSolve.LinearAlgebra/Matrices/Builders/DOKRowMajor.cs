using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: sorting the Dictionaries with Linq seems to take a lot of time.
namespace ISAAR.MSolve.LinearAlgebra.Matrices.Builders
{
    /// <summary>
    /// Use this class for building a sparse matrix, e.g. <see cref="CSCMatrix"/> not for operations. Convert to other matrix 
    /// formats once finished and use them instead for matrix operations. Only the non zero entries of the upper triangle are  
    /// stored. This class is optimized for building global positive definite matrices, where there is at least 1 entry per 
    /// column.
    /// </summary>
    public class DOKRowMajor: ISparseMatrix, IMatrixBuilder
    {
        /// <summary>
        /// See the rant in <see cref="DOKSymmetricColMajor.columns"/> about performance.
        /// </summary>
        private readonly Dictionary<int, double>[] rows;

        private DOKRowMajor(int numRows, int numCols, Dictionary<int, double>[] rows)
        {
            this.rows = rows;
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
                if (rows[rowIdx].TryGetValue(colIdx, out double val)) return val;
                else return 0.0;
            }
            set //not thread safe
            {
                rows[rowIdx][colIdx] = value;
            }
        }

        public static DOKRowMajor CreateEmpty(int numRows, int numCols)
        {
            var rows = new Dictionary<int, double>[numRows];
            for (int i = 0; i < numRows; ++i) rows[i] = new Dictionary<int, double>(); //Initial capacity may be optimized.
            return new DOKRowMajor(numRows, numCols, rows);
        }

        public static DOKRowMajor CreateIdentity(int order)
        {
            var rows = new Dictionary<int, double>[order];
            for (int j = 0; j < order; ++j)
            {
                var idenityRow = new Dictionary<int, double>(); //Initial capacity may be optimized.
                idenityRow[j] = 1.0;
                rows[j] = idenityRow;
            }
            return new DOKRowMajor(order, order, rows);
        }

        public static DOKRowMajor CreateFromSparseMatrix(ISparseMatrix matrix)
        {
            return CreateFromSparsePattern(matrix.NumRows, matrix.NumColumns, matrix.EnumerateNonZeros());
        }

        public static DOKRowMajor CreateFromSparsePattern(int numRows, int numColumns,
            IEnumerable<(int row, int col, double value)> nonZeroEntries)
        {
            DOKRowMajor dok = CreateEmpty(numRows, numColumns);
            foreach (var (row, col, val) in nonZeroEntries) dok.rows[row].Add(col, val);
            return dok;
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
            if (rows[rowIdx].TryGetValue(colIdx, out double oldValue))
            {
                rows[rowIdx][colIdx] = value + oldValue;
            }
            else rows[rowIdx][colIdx] = value;
            //The Dictionary rowss[rowIdx] is indexed twice in both cases. Is it possible to only index it once?
        }

        public void AddSubmatrix(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols)
        {
            foreach (var rowPair in subRowsToGlobalRows)
            {
                int subRow = rowPair.Key;
                int globalRow = rowPair.Value;
                foreach (var colPair in subColsToGlobalCols)
                {
                    int subCol = colPair.Key;
                    int globalCol = colPair.Value;
                    double elementValue = subMatrix[subRow, subCol];
                    if (rows[globalRow].TryGetValue(globalCol, out double oldGlobalValue))
                    {
                        rows[globalRow][globalCol] = elementValue + oldGlobalValue;
                    }
                    else rows[globalRow][globalCol] = elementValue;
                }
            }
        }

        /// <summary>
        /// Use this method if 1) both the global DOK and the element matrix are symmetric and 2) the rows and columns correspond
        /// to the same degrees of freedom. The caller is responsible for making sure that both matrices are symmetric and that 
        /// the dimensions and dofs of the element matrix and dof mappings match.
        /// </summary>
        /// <param name="subMatrix">The element submatrix, entries of which will be added to the global DOK. It must be 
        ///     symmetric and its <see cref="IIndexable2D.NumColumns"/> = <see cref="IIndexable2D.NumRows"/> must be equal to
        ///     elemenDofs.Length = globalDofs.Length.</param>
        /// <param name="subDofs">The entries in the element matrix to be added to the global matrix. Specificaly, pairs of 
        ///     (elementDofs[i], elementDofs[j]) will be added to (globalDofs[i], globalDofs[j]).</param>
        /// <param name="globalDofs">The entries in the global matrix where element matrix entries will be added to. Specificaly,
        ///     pairs of (elementDofs[i], elementDofs[j]) will be added to (globalDofs[i], globalDofs[j]).</param>
        public void AddSubmatrixSymmetric(IIndexable2D subMatrix, int[] subDofs, int[] globalDofs)
        {
            int n = subDofs.Length;
            for (int i = 0; i < n; ++i)
            {
                int subRow = subDofs[i];
                int globalRow = globalDofs[i];

                //Diagonal entry
                if (rows[globalRow].TryGetValue(globalRow, out double oldGlobalDiagValue))
                {
                    rows[globalRow][globalRow] = subMatrix[subRow, subRow] + oldGlobalDiagValue;
                }
                else rows[globalRow][globalRow] = subMatrix[subRow, subRow];

                //Non diagonal entries
                for (int j = 0; j < i; ++j)
                {
                    int subCol = subDofs[j];
                    int globalCol = globalDofs[j];
                    double newGlobalValue = subMatrix[subRow, subRow];
                    // Only check the lower triangle. If the DOK matrix is not symmetric, this will cause errors
                    if (rows[globalRow].TryGetValue(globalCol, out double oldGlobalValue))
                    {
                        newGlobalValue += oldGlobalValue;
                    }
                    rows[globalRow][globalCol] = newGlobalValue;
                    rows[globalCol][globalRow] = newGlobalValue;
                }
            }
        }
        #endregion

        /// <summary>
        /// Creates the CSC arrays.
        /// </summary>
        /// <param name="sortColsOfEachCol">True to sort the column indices of the CSR matrix between rowOffsets[i] and 
        ///     rowOffsets[i+1] in ascending order. False to leave them unordered. Ordered columns might result in better  
        ///     performance during multiplications or they might be required by 3rd party libraries. Leaving them unordered will 
        ///     be faster during creation of the CSR matrix.</param>
        /// <returns></returns>
        public (double[] values, int[] rowIndices, int[] columnOffsets) BuildCSRArrays(bool sortColsOfEachCol)
        {
            int[] rowOffsets = new int[NumRows + 1];
            int nnz = 0;
            for (int i = 0; i < NumRows; ++i)
            {
                rowOffsets[i] = nnz;
                nnz += rows[i].Count;
            }
            rowOffsets[NumRows] = nnz; //The last CSR entry is nnz.

            int[] colIndices = new int[nnz];
            double[] values = new double[nnz];
            int counter = 0;
            for (int i = 0; i < NumRows; ++i)
            {
                if (sortColsOfEachCol)
                {
                    foreach (var rowVal in rows[i].OrderBy(pair => pair.Key))
                    {
                        colIndices[counter] = rowVal.Key;
                        values[counter] = rowVal.Value;
                        ++counter;
                    }
                }
                else
                {
                    foreach (var rowVal in rows[i])
                    {
                        colIndices[counter] = rowVal.Key;
                        values[counter] = rowVal.Value;
                        ++counter;
                    }
                }
            }

            return (values, colIndices, rowOffsets);
        }

        /// <summary>
        /// Creates a CSC matrix.
        /// </summary>
        /// <param name="sortColsOfEachCol">True to sort the column indices of the CSR matrix between rowOffsets[i] and 
        ///     rowOffsets[i+1] in ascending order. False to leave them unordered. Ordered columns might result in better  
        ///     performance during multiplications or they might be required by 3rd party libraries. Leaving them unordered will 
        ///     be faster during creation of the CSR matrix.</param>
        /// <returns></returns>
        public CSRMatrix BuildCSRMatrix(bool sortColsOfEachCol)
        {
            (double[] values, int[] colIndices, int[] rowOffsets) = BuildCSRArrays(sortColsOfEachCol);
            return CSRMatrix.CreateFromArrays(NumRows, NumColumns, values, colIndices, rowOffsets, false);
        }

        public int CountNonZeros()
        {
            int count = 0;
            for (int j = 0; j < NumColumns; ++j) count += rows[j].Count;
            return count;
        }

        public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
        {
            for (int i = 0; i < NumRows; ++i)
            {
                foreach (var colVal in rows[i])
                {
                    yield return (i, colVal.Key, colVal.Value);
                }
            }
        }

        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
        {
            return DenseStrategies.AreEqual(this, other, tolerance);
        }

        /// <summary>
        /// Returns the diagonal as <see cref="Vector"/> and the index of the first zero entry. If there are no zero entries,
        /// -1 is returned as the index.
        /// </summary>
        /// <returns></returns>
        public (Vector diagonal, int firstZeroIdx) GetDiagonal()
        {
            (double[] diagonal, int firstZeroIdx) = GetDiagonalAsArray();
            return (Vector.CreateFromArray(diagonal, false), firstZeroIdx);
        }

        /// <summary>
        /// Returns the diagonal as <see cref="double[]"/> and the index of the first zero entry. If there are no zero entries,
        /// -1 is returned as the index.
        /// </summary>
        /// <returns></returns>
        public (double[] diagonal, int firstZeroIdx) GetDiagonalAsArray()
        {
            Preconditions.CheckSquare(this);
            double[] diag = new double[NumRows];
            int firstZeroIdx = -1;
            for (int i = 0; i < NumRows; ++i)
            {
                bool isStored = rows[i].TryGetValue(i, out double val);
                if (isStored) diag[i] = val;
                else
                {
                    diag[i] = 0.0;
                    firstZeroIdx = i;
                }
            }
            return (diag, firstZeroIdx);
        }

        public SparseFormat GetSparseFormat()
        {
            (double[] values, int[] colIndices, int[] rowOffsets) = BuildCSRArrays(false);
            var format = new SparseFormat();
            format.RawValuesTitle = "Values";
            format.RawValuesArray = values;
            format.RawIndexArrays.Add("Column indices", colIndices);
            format.RawIndexArrays.Add("Row offsets", rowOffsets);
            return format;
        }

        /// <summary>
        /// At least as fast as creating a CSR matrix and multiplying it with 1 vector. If more than 1 multiplications are 
        /// needed, then the CSR matrix should be created.
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="avoidBuilding">If true, no matrices will be built internally. If false, such matrices may be built if 
        ///     the method decides that using them is faster (e.g. fast native library), however that requires memory for 
        ///     building them, which may not be available.</param>
        /// <returns></returns>
        public Vector MultiplyRight(Vector vector, bool avoidBuilding = false)
        {
            // MKL >> managed code. Just don't sort the CSR.
            bool buildCSR = (!avoidBuilding) && CSRMatrix.UseMKL; 
            if (buildCSR) return BuildCSRMatrix(false).MultiplyRight(vector); 

            Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
            var result = new double[this.NumRows];
            for (int i = 0; i < NumRows; ++i)
            {
                double dot = 0.0;
                foreach (var colValPair in rows[i]) dot += colValPair.Value * vector[colValPair.Key];
                result[i] = dot;
            }
            return Vector.CreateFromArray(result, false);
        }
    }
}
