using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices.Builders
{
    /// <summary>
    /// Use this class for building large sparse matrices, e.g. <see cref="CscMatrix"/>, not for operations. Convert to other 
    /// matrix formats once finished and use them instead for matrix operations. The large matrices and their properties will be 
    /// characterized as "global" in this namespace. This class is optimized for building global matrices with at least 1 entry 
    /// per column and column major storage formats.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class DokColMajor: ISparseMatrix, IGeneralMatrixBuilder
    {
        /// <summary>
        /// See the rant in <see cref="DokSymmetric.columns"/> about performance.
        /// </summary>
        private Dictionary<int, double>[] columns;

        private DokColMajor(int numRows, int numCols, Dictionary<int, double>[] columns)
        {
            this.columns = columns;
            this.NumRows = numRows;
            this.NumColumns = numCols;
        }

        /// <summary>
        /// The number of columns of the matrix.
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get; }

        /// <summary>
        /// See <see cref="IIndexable2D.this[int, int]"/>, <see cref="IMatrixBuilder.this[int, int]"/>.
        /// </summary>
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

        /// <summary>
        /// Initializes a new instance of <see cref="DokColMajor"/> with the specified matrix dimensions and all entries being
        /// equal to 0.
        /// </summary>
        /// <param name="numRows">The number of rows of the matrix to build.</param>
        /// <param name="numCols">The number of columns of the matrix to build.</param>
        public static DokColMajor CreateEmpty(int numRows, int numCols)
        {
            var columns = new Dictionary<int, double>[numCols];
            for (int j = 0; j < numCols; ++j) columns[j] = new Dictionary<int, double>(); //Initial capacity may be optimized.
            return new DokColMajor(numRows, numCols, columns);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="DokColMajor"/> with the specified matrix dimensions and entries being 
        /// the same as the identity matrix.
        /// </summary>
        /// <param name="order">The number of rows/columns of the square matrix to build.</param>
        public static DokColMajor CreateIdentity(int order)
        {
            var columns = new Dictionary<int, double>[order];
            for (int j = 0; j < order; ++j)
            {
                var idenityCol = new Dictionary<int, double>(); //Initial capacity may be optimized.
                idenityCol[j] = 1.0;
                columns[j] = idenityCol;
            }
            return new DokColMajor(order, order, columns);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="DokColMajor"/> with the specified matrix dimensions and the non-zero 
        /// entries of the provided sparse matrix.
        /// </summary>
        /// <param name="matrix">A sparse matrix whose dimensions and non-zero entries will be used to intialize the new
        ///     <see cref="DokColMajor"/>.</param>
        public static DokColMajor CreateFromSparseMatrix(ISparseMatrix matrix)
            => CreateFromSparsePattern(matrix.NumRows, matrix.NumColumns, matrix.EnumerateNonZeros());

        /// <summary>
        /// Initializes a new instance of <see cref="DokColMajor"/> with the specified matrix dimensions and the non-zero 
        /// entries defined by the provided pattern.
        /// </summary>
        /// <param name="numRows">The number of rows of the matrix to build.</param>
        /// <param name="numCols">The number of columns of the matrix to build.</param>
        /// <param name="nonZeroEntries">The non-zero entries of the matrix to build.</param>
        public static DokColMajor CreateFromSparsePattern(int numRows, int numColumns, 
            IEnumerable<(int row, int col, double value)> nonZeroEntries)
        {
            DokColMajor dok = CreateEmpty(numRows, numColumns);
            foreach (var (row, col, val) in nonZeroEntries) dok.columns[col].Add(row, val);
            return dok;
        }

        /// <summary>
        /// See <see cref="IMatrixBuilder.AddToEntry(int, int, double)"/>
        /// </summary>
        public void AddToEntry(int rowIdx, int colIdx, double value)
        {
            if (columns[colIdx].TryGetValue(rowIdx, out double oldValue))
            {
                columns[colIdx][rowIdx] = value + oldValue;
            }
            else columns[colIdx][rowIdx] = value;
            //The Dictionary columns[rowIdx] is indexed twice in both cases. Is it possible to only index it once?
        }

        /// <summary>
        /// See <see cref="IGeneralMatrixBuilder.AddSubmatrix(IIndexable2D, IReadOnlyDictionary{int, int}, 
        /// IReadOnlyDictionary{int, int})"/>.
        /// </summary>
        public void AddSubmatrix(IIndexable2D subMatrix, 
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols)
        {
            foreach (var colPair in subColsToGlobalCols) // Col major ordering is the default
            {
                int subCol = colPair.Key;
                int globalCol = colPair.Value;
                foreach (var rowPair in subRowsToGlobalRows)
                {
                    int subRow = rowPair.Key;
                    int globalRow = rowPair.Value;
                    double subValue = subMatrix[subRow, subCol];
                    if (columns[globalCol].TryGetValue(globalRow, out double oldGlobalValue))
                    {
                        columns[globalCol][globalRow] = subValue + oldGlobalValue;
                    }
                    else columns[globalCol][globalRow] = subValue;
                }
            }
        }

        /// <summary>
        /// See <see cref="IGeneralMatrixBuilder.AddSubmatrix(IIndexable2D, int[], int[], int[], int[])"/>.
        /// </summary>
        public void AddSubmatrix(IIndexable2D subMatrix, int[] subMatrixRows, int[] globalMatrixRows, 
            int[] subMatrixCols, int[] globalMatrixCols)
        {
            int numRows = subMatrixRows.Length;
            int numCols = subMatrixCols.Length;
            Debug.Assert(numRows == globalMatrixRows.Length);
            Debug.Assert(numCols == globalMatrixCols.Length);

            for (int j = 0; j < numCols; ++j)
            {
                int subCol = subMatrixCols[j];
                int globalCol = globalMatrixCols[j];
                Debug.Assert((globalCol >= 0) && (globalCol < NumColumns));

                for (int i = 0; i < numRows; ++i)
                {
                    int subRow = subMatrixRows[i];
                    int globalRow = globalMatrixRows[i];
                    Debug.Assert((globalRow >= 0) && (globalRow < NumRows));

                    double subVal = subMatrix[subRow, subCol];
                    columns[globalCol].TryGetValue(globalRow, out double oldGlobalVal); 
                    columns[globalCol][globalRow] = subVal + oldGlobalVal; // default old value = 0.0, if the entry is new
                }
            }
        }

        /// <summary>
        /// See <see cref="IGeneralMatrixBuilder.AddSubmatrixSymmetric(IIndexable2D, IReadOnlyDictionary{int, int})"/>.
        /// </summary>
        public void AddSubmatrixSymmetric(IIndexable2D subMatrix, IReadOnlyDictionary<int, int> subIndicesToGlobalIndices)
        {
            foreach (var colPair in subIndicesToGlobalIndices)
            {
                int subCol = colPair.Key;
                int globalCol = colPair.Value;
                foreach (var rowPair in subIndicesToGlobalIndices)
                {
                    int subRow = rowPair.Key;
                    int globalRow = rowPair.Value;

                    if (globalCol > globalRow)
                    {
                        double subValue = subMatrix[subRow, subCol];
                        columns[globalCol].TryGetValue(globalRow, out double oldGlobalValue); // default value = 0.0
                        double newGlobalValue = oldGlobalValue + subValue;
                        columns[globalCol][globalRow] = newGlobalValue;
                        columns[globalRow][globalCol] = newGlobalValue;
                    }
                    else if (globalCol == globalRow)
                    {
                        double subValue = subMatrix[subRow, subCol];
                        columns[globalCol].TryGetValue(globalRow, out double oldGlobalValue); // default value = 0.0
                        columns[globalCol][globalRow] = oldGlobalValue + subValue;
                    }
                }
            }
        }

        /// <summary>
        /// See <see cref="IGeneralMatrixBuilder.AddSubmatrixSymmetric(IIndexable2D, int[], int[])"/>.
        /// </summary>
        public void AddSubmatrixSymmetric(IIndexable2D subMatrix, int[] subMatrixIndices, int[] globalIndices)
        {
            Debug.Assert(subMatrix.NumRows == subMatrix.NumColumns);
            Debug.Assert(globalIndices.Length == subMatrixIndices.Length);

            int numRelevantRows = subMatrixIndices.Length;
            for (int j = 0; j < numRelevantRows; ++j)
            {
                int subCol = subMatrixIndices[j];
                int globalCol = globalIndices[j];
                for (int i = 0; i < numRelevantRows; ++i)
                {
                    int subRow = subMatrixIndices[i];
                    int globalRow = globalIndices[i];

                    if (globalCol > globalRow)
                    {
                        double subValue = subMatrix[subRow, subCol];
                        columns[globalCol].TryGetValue(globalRow, out double oldGlobalValue); // default value = 0.0
                        double newGlobalValue = oldGlobalValue + subValue;
                        columns[globalCol][globalRow] = newGlobalValue;
                        columns[globalRow][globalCol] = newGlobalValue;
                    }
                    else if (globalCol == globalRow)
                    {
                        double subValue = subMatrix[subRow, subCol];
                        columns[globalCol].TryGetValue(globalRow, out double oldGlobalValue); // default value = 0.0
                        columns[globalCol][globalRow] = oldGlobalValue + subValue;
                    }
                }
            }
        }

        /// <summary>
        /// Creates the values and indexing arrays in CSC storage format of the current matrix. This method should be 
        /// called after fully defining the matrix in <see cref="DokColMajor"/> format.
        /// </summary>
        /// <param name="sortColsOfEachCol">True to sort the column indices of the CSC matrix between colOffsets[j] and 
        ///     colOffsets[j+1] in ascending order. False to leave them unordered. Ordered rows might result in better  
        ///     performance during multiplications or they might be required by 3rd party libraries. Conversely, leaving them 
        ///     unordered will be faster during creation of the CSC matrix.</param>
        /// <exception cref="EmptyMatrixBuilderException">Thrown if no non-zero entries have been defined yet.</exception>
        public (double[] values, int[] rowIndices, int[] colOffsets) BuildCscArrays(bool sortRowsOfEachCol)
        {
            int[] colOffsets = new int[NumColumns + 1];
            int nnz = 0;
            for (int j = 0; j < NumColumns; ++j)
            {
                colOffsets[j] = nnz;
                nnz += columns[j].Count;
            }
            if (nnz == 0) throw new EmptyMatrixBuilderException("Cannot build CSC arrays from a DOK with nnz = 0.");
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
        /// Initializes a <see cref="CscMatrix"/> representation of the current matrix. This method should be 
        /// called after fully defining the matrix in <see cref="DokColMajor"/> format.
        /// </summary>
        /// <param name="sortColsOfEachCol">True to sort the column indices of the CSC matrix between colOffsets[j] and 
        ///     colOffsets[j+1] in ascending order. False to leave them unordered. Ordered rows might result in better  
        ///     performance during multiplications or they might be required by 3rd party libraries. Conversely, leaving them 
        ///     unordered will be faster during creation of the CSC matrix.</param>
        /// <exception cref="EmptyMatrixBuilderException">Thrown if no non-zero entries have been defined yet.</exception>
        public CscMatrix BuildCscMatrix(bool sortRowsOfEachCol)
        {
            (double[] values, int[] rowIndices, int[] colOffsets) = BuildCscArrays(sortRowsOfEachCol);
            return CscMatrix.CreateFromArrays(NumRows, NumColumns, values, rowIndices, colOffsets, false);
        }

        /// <summary>
        /// Frees all memory held by this <see cref="DokColMajor"/> instance. Afterwards this object cannot be reused. Therefore 
        /// this method should be called to save up space after building the matrix or its internal arrays.
        /// </summary>
        public void Clear() => columns = null;

        /// <summary>
        /// See <see cref="ISparseMatrix.CountNonZeros"/>.
        /// </summary>
        public int CountNonZeros()
        {
            int count = 0;
            for (int j = 0; j < NumColumns; ++j) count += columns[j].Count;
            return count;
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.EnumerateNonZeros"/>.
        /// </summary>
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

        /// <summary>
        /// See <see cref="IIndexable2D.Equals(IIndexable2D, double)"/>.
        /// </summary>
        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
        {
            return DenseStrategies.AreEqual(this, other, tolerance);
        }

        /// <summary>
        /// Returns a <see cref="Vector"/> with the entries of the matrix's main diagonal and the index of the first zero entry.
        /// If there are no zero entries, -1 is returned as the index. The <see cref="DokColMajor"/> matrix must be square.
        /// </summary>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
        public (Vector diagonal, int firstZeroIdx) GetDiagonal()
        {
            (double[] diagonal, int firstZeroIdx) = GetDiagonalAsArray();
            return (Vector.CreateFromArray(diagonal, false), firstZeroIdx);
        }

        /// <summary>
        /// Returns an array with the entries of the matrix's main diagonal and the index of the first zero entry.
        /// If there are no zero entries, -1 is returned as the index. The <see cref="DokColMajor"/> matrix must be square.
        /// </summary>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
        public (double[] diagonal, int firstZeroIdx) GetDiagonalAsArray()
        {
            Preconditions.CheckSquare(this);
            double[] diag = new double[NumColumns];
            int firstZeroIdx = -1;
            for (int j = 0; j < NumColumns; ++j)
            {
                bool isStored = columns[j].TryGetValue(j, out double val);
                if (isStored) diag[j] = val;
                else
                {
                    diag[j] = 0.0;
                    firstZeroIdx = j;
                }
            }
            return (diag, firstZeroIdx);
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.GetSparseFormat"/>.
        /// </summary>
        public SparseFormat GetSparseFormat()
        {
            (double[] values, int[] rowIndices, int[] colOffsets) = BuildCscArrays(false);
            var format = new SparseFormat();
            format.RawValuesTitle = "Values";
            format.RawValuesArray = values;
            format.RawIndexArrays.Add("Row indices", rowIndices);
            format.RawIndexArrays.Add("Column offsets", colOffsets);
            return format;
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
        private void AddSubmatrixSymmetricOLD(IIndexable2D subMatrix, int[] subDofs, int[] globalDofs) //TODO: this should be reworked
        {
            int n = subDofs.Length;
            for (int j = 0; j < n; ++j)
            {
                int subCol = subDofs[j];
                int globalCol = globalDofs[j];

                //Diagonal entry
                if (columns[globalCol].TryGetValue(globalCol, out double oldGlobalDiagValue))
                {
                    columns[globalCol][globalCol] = subMatrix[subCol, subCol] + oldGlobalDiagValue;
                }
                else columns[globalCol][globalCol] = subMatrix[subCol, subCol];

                //Non diagonal entries
                for (int i = 0; i < j; ++i)
                {
                    int subRow = subDofs[j];
                    int globalRow = globalDofs[j];
                    double newGlobalValue = subMatrix[subRow, subCol];
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
    }
}
