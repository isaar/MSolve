using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: sorting the Dictionaries with Linq seems to take a lot of time.
namespace ISAAR.MSolve.LinearAlgebra.Matrices.Builders
{
    /// <summary>
    /// Use this class for building large sparse matrices, e.g. <see cref="CsrMatrix"/>, not for operations. Convert to other 
    /// matrix formats once finished and use them instead for matrix operations. The large matrices and their properties will be 
    /// characterized as "global" in this namespace. This class is optimized for building global matrices with at least 1 entry 
    /// per column and row major storage formats.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class DokRowMajor: ISparseMatrix, IGeneralMatrixBuilder
    {
        /// <summary>
        /// See the rant in <see cref="DokSymmetric.columns"/> about performance.
        /// </summary>
        private Dictionary<int, double>[] rows;

        private DokRowMajor(int numRows, int numCols, Dictionary<int, double>[] rows)
        {
            this.rows = rows;
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
                if (rows[rowIdx].TryGetValue(colIdx, out double val)) return val;
                else return 0.0;
            }
            set //not thread safe
            {
                rows[rowIdx][colIdx] = value;
            }
        }

        /// <summary>
        /// Initializes a new instance of <see cref="DokRowMajor"/> with the specified matrix dimensions and all entries being
        /// equal to 0.
        /// </summary>
        /// <param name="numRows">The number of rows of the matrix to build.</param>
        /// <param name="numCols">The number of columns of the matrix to build.</param>
        public static DokRowMajor CreateEmpty(int numRows, int numCols)
        {
            var rows = new Dictionary<int, double>[numRows];
            for (int i = 0; i < numRows; ++i) rows[i] = new Dictionary<int, double>(); //Initial capacity may be optimized.
            return new DokRowMajor(numRows, numCols, rows);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="DokRowMajor"/> with the specified matrix dimensions and entries being 
        /// the same as the identity matrix.
        /// </summary>
        /// <param name="order">The number of rows/columns of the square matrix to build.</param>
        public static DokRowMajor CreateIdentity(int order)
        {
            var rows = new Dictionary<int, double>[order];
            for (int j = 0; j < order; ++j)
            {
                var idenityRow = new Dictionary<int, double>(); //Initial capacity may be optimized.
                idenityRow[j] = 1.0;
                rows[j] = idenityRow;
            }
            return new DokRowMajor(order, order, rows);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="DokRowMajor"/> with the specified matrix dimensions and the non-zero 
        /// entries of the provided sparse matrix.
        /// </summary>
        /// <param name="matrix">A sparse matrix whose dimensions and non-zero entries will be used to intialize the new
        ///     <see cref="DokRowMajor"/>.</param>
        public static DokRowMajor CreateFromSparseMatrix(ISparseMatrix matrix)
            => CreateFromSparsePattern(matrix.NumRows, matrix.NumColumns, matrix.EnumerateNonZeros());

        /// <summary>
        /// Initializes a new instance of <see cref="DokRowMajor"/> with the specified matrix dimensions and the non-zero 
        /// entries defined by the provided pattern.
        /// </summary>
        /// <param name="numRows">The number of rows of the matrix to build.</param>
        /// <param name="numCols">The number of columns of the matrix to build.</param>
        /// <param name="nonZeroEntries">The non-zero entries of the matrix to build.</param>
        public static DokRowMajor CreateFromSparsePattern(int numRows, int numColumns,
            IEnumerable<(int row, int col, double value)> nonZeroEntries)
        {
            DokRowMajor dok = CreateEmpty(numRows, numColumns);
            foreach (var (row, col, val) in nonZeroEntries) dok.rows[row].Add(col, val);
            return dok;
        }

        /// <summary>
        /// See <see cref="IMatrixBuilder.AddToEntry(int, int, double)"/>
        /// </summary>
        public void AddToEntry(int rowIdx, int colIdx, double value)
        {
            if (rows[rowIdx].TryGetValue(colIdx, out double oldValue))
            {
                rows[rowIdx][colIdx] = value + oldValue;
            }
            else rows[rowIdx][colIdx] = value;
            //The Dictionary rowss[rowIdx] is indexed twice in both cases. Is it possible to only index it once?
        }

        /// <summary>
        /// See <see cref="IGeneralMatrixBuilder.AddSubmatrix(IIndexable2D, IReadOnlyDictionary{int, int}, 
        /// IReadOnlyDictionary{int, int})"/>.
        /// </summary>
        public void AddSubmatrix(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols)
        {
            foreach (var rowPair in subRowsToGlobalRows)
            {
                int subRow = rowPair.Key;
                int globalRow = rowPair.Value;
                //Debug.Assert((globalRow >= 0) && (globalRow < NumRows));
                foreach (var colPair in subColsToGlobalCols)
                {
                    int subCol = colPair.Key;
                    int globalCol = colPair.Value;
                    //Debug.Assert((globalCol >= 0) && (globalCol < NumColumns));
                    double subValue = subMatrix[subRow, subCol];
                    if (rows[globalRow].TryGetValue(globalCol, out double oldGlobalValue))
                    {
                        rows[globalRow][globalCol] = subValue + oldGlobalValue;
                    }
                    else rows[globalRow][globalCol] = subValue;
                }
            }
        }

        /// <summary>
        /// See <see cref="IGeneralMatrixBuilder.AddSubmatrixSymmetric(IIndexable2D, IReadOnlyDictionary{int, int})"/>.
        /// </summary>
        public void AddSubmatrixSymmetric(IIndexable2D subMatrix, IReadOnlyDictionary<int, int> subIndicesToGlobalIndices)
        {
            foreach (var rowPair in subIndicesToGlobalIndices)
            {
                int subRow = rowPair.Key;
                int globalRow = rowPair.Value;
                foreach (var colPair in subIndicesToGlobalIndices)
                {
                    int subCol = colPair.Key;
                    int globalCol = colPair.Value;

                    if (globalCol > globalRow)
                    {
                        double subValue = subMatrix[subRow, subCol];
                        rows[globalRow].TryGetValue(globalCol, out double oldGlobalValue); // default value = 0.0
                        double newGlobalValue = oldGlobalValue + subValue;
                        rows[globalRow][globalCol] = newGlobalValue;
                        rows[globalCol][globalRow] = newGlobalValue;
                    }
                    else if (globalCol == globalRow)
                    {
                        double subValue = subMatrix[subRow, subCol];
                        rows[globalRow].TryGetValue(globalCol, out double oldGlobalValue); // default value = 0.0
                        rows[globalRow][globalCol] = oldGlobalValue + subValue;
                    }
                }
            }
        }

        /// <summary>
        /// Creates the values and indexing arrays in CSR storage format of the current matrix. This method should be 
        /// called after fully defining the matrix in <see cref="DokRowMajor"/> format.
        /// </summary>
        /// <param name="sortColsOfEachCol">True to sort the column indices of the CSR matrix between rowOffsets[i] and 
        ///     rowOffsets[i+1] in ascending order. False to leave them unordered. Ordered columns might result in better  
        ///     performance during multiplications or they might be required by 3rd party libraries. Conversely, leaving them 
        ///     unordered will be faster during creation of the CSR matrix.</param>
        /// <exception cref="EmptyMatrixBuilderException">Thrown if no non-zero entries have been defined yet.</exception>
        public (double[] values, int[] colIndices, int[] rowOffsets) BuildCsrArrays(bool sortColsOfEachCol)
        {
            int[] rowOffsets = new int[NumRows + 1];
            int nnz = 0;
            for (int i = 0; i < NumRows; ++i)
            {
                rowOffsets[i] = nnz;
                nnz += rows[i].Count;
            }
            if (nnz == 0) throw new EmptyMatrixBuilderException("Cannot build CSR arrays from a DOK with nnz = 0.");
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
        /// Initializes a <see cref="CsrMatrix"/> representation of the current matrix. This method should be 
        /// called after fully defining the matrix in <see cref="DokRowMajor"/> format.
        /// </summary>
        /// <param name="sortColsOfEachRow">True to sort the column indices of the CSR matrix between rowOffsets[i] and 
        ///     rowOffsets[i+1] in ascending order. False to leave them unordered. Ordered columns might result in better  
        ///     performance during multiplications or they might be required by 3rd party libraries. Conversely, leaving them 
        ///     unordered will be faster during creation of the CSR matrix.</param>
        /// <exception cref="EmptyMatrixBuilderException">Thrown if no non-zero entries have been defined yet.</exception>
        public CsrMatrix BuildCsrMatrix(bool sortColsOfEachRow)
        {
            (double[] values, int[] colIndices, int[] rowOffsets) = BuildCsrArrays(sortColsOfEachRow);
            return CsrMatrix.CreateFromArrays(NumRows, NumColumns, values, colIndices, rowOffsets, false);
        }

        /// <summary>
        /// Frees all memory held by this <see cref="DokRowMajor"/> instance. Afterwards this object cannot be reused. Therefore 
        /// this method should be called to save up space after building the matrix or its internal arrays.
        /// </summary>
        public void Clear() => rows = null;

        /// <summary>
        /// See <see cref="ISparseMatrix.CountNonZeros"/>.
        /// </summary>
        public int CountNonZeros()
        {
            int count = 0;
            for (int j = 0; j < NumColumns; ++j) count += rows[j].Count;
            return count;
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.EnumerateNonZeros"/>.
        /// </summary>
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

        /// <summary>
        /// See <see cref="IIndexable2D.Equals(IIndexable2D, double)"/>.
        /// </summary>
        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
        {
            return DenseStrategies.AreEqual(this, other, tolerance);
        }

        /// <summary>
        /// Returns a <see cref="Vector"/> with the entries of the matrix's main diagonal and the index of the first zero entry.
        /// If there are no zero entries, -1 is returned as the index. The <see cref="DokRowMajor"/> matrix must be square.
        /// </summary>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
        public (Vector diagonal, int firstZeroIdx) GetDiagonal()
        {
            (double[] diagonal, int firstZeroIdx) = GetDiagonalAsArray();
            return (Vector.CreateFromArray(diagonal, false), firstZeroIdx);
        }

        /// <summary>
        /// Returns an array with the entries of the matrix's main diagonal and the index of the first zero entry.
        /// If there are no zero entries, -1 is returned as the index. The <see cref="DokRowMajor"/> matrix must be square.
        /// </summary>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
        public (double[] diagonal, int firstZeroIdx) GetDiagonalAsArray() //TODO: the -1 sentinel value should be a constant somewhere
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

        /// <summary>
        /// See <see cref="ISparseMatrix.GetSparseFormat"/>.
        /// </summary>
        public SparseFormat GetSparseFormat()
        {
            (double[] values, int[] colIndices, int[] rowOffsets) = BuildCsrArrays(false);
            var format = new SparseFormat();
            format.RawValuesTitle = "Values";
            format.RawValuesArray = values;
            format.RawIndexArrays.Add("Column indices", colIndices);
            format.RawIndexArrays.Add("Row offsets", rowOffsets);
            return format;
        }

        /// <summary>
        /// Performs the matrix-vector multiplication: this * <paramref name="vector"/>. This method is at least as fast as 
        /// creating a CSR matrix and multiplying it with 1 vector. If more than 1 multiplications are needed, then the CSR 
        /// matrix should be created.
        /// </summary>
        /// <param name="vector">A vector with <see cref="IIndexable1D.Length"/> being equal to this.<see cref="NumColumns"/>.
        ///     </param>
        /// <param name="avoidBuilding">If true, no matrices will be built internally. If false, such matrices may be built if 
        ///     the method decides that using them is faster (e.g. fast native library), however that requires memory for 
        ///     building them, which may not be available.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if 
        ///     <paramref name="vector"/>.<see cref="IIndexable1D.Length"/> != this.<see cref="NumColumns"/>.</exception>
        public Vector MultiplyRight(Vector vector, bool avoidBuilding = false)
        {
            // MKL functions are way faster than managed code. Just don't sort the CSR.
            bool buildCSR = (!avoidBuilding) && CsrMatrix.UseMKL; 
            if (buildCSR) return BuildCsrMatrix(false).MultiplyRight(vector); 

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
        private void AddSubmatrixSymmetric(IIndexable2D subMatrix, int[] subDofs, int[] globalDofs) //TODO: this should be reworked
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
    }
}
