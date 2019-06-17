using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.CompilerServices;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Add a create from dense method to facilitate testing.
//TODO: Also provide a method GetSubmatrixWithout that would remove rows/cols. Actually that would be better as Remove(int[])
//      in a DOK with a Dictionary<Dictionary<int, double>> as a backing storage to easily remove stuff.
namespace ISAAR.MSolve.LinearAlgebra.Matrices.Builders
{
    /// <summary>
    /// Use this class for building large symmetric sparse matrices not for operations. Convert to other matrix formats once 
    /// finished and use them instead for matrix operations. The large matrices and their properties will be characterized as 
    /// "global" in this namespace. This class is optimized for building symmetric global matrices, where only entries of the
    /// upper triangle are explicitly stored, the storage format is column major and there is at least 1 entry per column.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class DokSymmetric: ISparseSymmetricMatrix, ISymmetricMatrixBuilder
    {
        /// <summary>
        /// An array of dictionaries is more efficent and perhaps easier to work with than a dictionary of dictionaries. There 
        /// is usually at least 1 non zero entry in each column. Otherwise this data structure wastes space, but if there were  
        /// many empty rows, perhaps another matrix format is more appropriate.
        /// To get the value: columns[colIdx][rowIdx] = value. 
        /// To get the row-value subdictionary: columns[colIdx] = Dictionary[int, double]
        /// TODO: Things to benchmark (keys = row indices): 
        /// 1) Dictionary + sort on converting: copy unordered keys and values to arrays and sort both arrays according to keys 
        /// 2) Dictionary + sort on converting: Use LINQ to sort key-value pairs and then copy each one to the arrays
        /// 3) SortedDictionary: O(log(n)) insertion but already ordered. Some entries are inserted more than once though!
        /// However this is wasteful if the sparse matrix doesn't need to be ordered. Thus I should use IDictionary and provide 
        /// static methods for Skyline (descending keys), CSC (ascending keys) or unordered (can use simple Dictionary)
        /// Perhaps a Dictionary should be used instead of SortedDictionary and only sort each column independently before 
        /// converting.
        /// </summary>
        private Dictionary<int, double>[] columns;
        private readonly int order;

        private DokSymmetric(int order, Dictionary<int, double>[] columns)
        {
            this.columns = columns;
            this.order = order;
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
                if (columns[colIdx].TryGetValue(rowIdx, out double val)) return val;
                else return 0.0;
            }
            set //not thread safe
            {
                ProcessIndices(ref rowIdx, ref colIdx);
                columns[colIdx][rowIdx] = value;
            }
        }

        /// <summary>
        /// Initializes a new instance of <see cref="DokColMajor"/> with the specified matrix dimensions and all entries being
        /// equal to 0.
        /// </summary>
        /// <param name="order">The number of rows/columns of the symmetric matrix to build.</param>
        public static DokSymmetric CreateEmpty(int order)
        {
            var columns = new Dictionary<int, double>[order];
            for (int j = 0; j < order; ++j) columns[j] = new Dictionary<int, double>(); //Initial capacity may be optimized.
            return new DokSymmetric(order, columns);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="DokSymmetric"/> with the specified matrix dimensions and entries 
        /// being the same as the identity matrix.
        /// WARNING: if you initilize the DOK to identity, then you need to overwrite the diagonal entries with the new values 
        /// (in FEM you only write to each diagonal entry once), instead of adding to them. This precludes the use of the methods 
        /// that add a whole submatrix to the DOK. Therefore, it is more convenient to call <see cref="CreateEmpty(int)"/> 
        /// and after all submatrices have been added, call <see cref="SetStructuralZeroDiagonalEntriesToUnity"/>.
        /// </summary>
        /// <param name="order">The number of rows/columns of the symmetric matrix to build.</param>
        public static DokSymmetric CreateIdentity(int order)
        {
            var columns = new Dictionary<int, double>[order];
            for (int j = 0; j < order; ++j) 
            {
                var idenityCol = new Dictionary<int, double>(); //Initial capacity may be optimized.
                idenityCol[j] = 1.0;
                columns[j] = idenityCol; 
            }
            return new DokSymmetric(order, columns);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="DokSymmetric"/> with the specified matrix dimensions and the 
        /// non-zero entries of the provided sparse matrix.
        /// </summary>
        /// <param name="matrix">A symmetric sparse matrix whose dimensions and non-zero entries will be used to intialize the 
        ///     new <see cref="DokSymmetric"/>. Its symmetry will not be checked explicitly.</param>
        public static DokSymmetric CreateFromSparseMatrix(ISparseMatrix matrix)
        {
            Preconditions.CheckSquare(matrix);
            return CreateFromSparsePattern(matrix.NumColumns, matrix.EnumerateNonZeros());
        }

        /// <summary>
        /// Initializes a new instance of <see cref="DokSymmetric"/> with the specified matrix dimensions and the 
        /// non-zero entries defined by the provided pattern.
        /// </summary>
        /// <param name="numRows">The number of rows/columns of the symmetric matrix to build.</param>
        /// <param name="nonZeroEntries">The non-zero entries of the symmetric matrix to build. Its symmetry will not be 
        ///     checked explicitly.</param>
        public static DokSymmetric CreateFromSparsePattern(int order, 
            IEnumerable<(int row, int col, double value)> nonZeroEntries)
        {
            DokSymmetric dok = CreateEmpty(order);
            foreach (var (row, col, val) in nonZeroEntries)
            {
                if ( row <= col) dok.columns[col][row] = val;
            }
            return dok;
        }

        /// <summary>
        /// Initializes a new instance of <see cref="DokSymmetric"/> with the specified matrix dimensions and the 
        /// non-zero entries of the provided symmetric sparse matrix.
        /// </summary>
        /// <param name="matrix">A symmetric sparse matrix whose dimensions and non-zero entries will be used to intialize the 
        ///     new <see cref="DokSymmetric"/>.</param>
        public static DokSymmetric CreateFromSparseSymmetricMatrix(ISparseSymmetricMatrix matrix)
        {
            Preconditions.CheckSquare(matrix); //This should not be needed
            return CreateFromSparseSymmetricPattern(matrix.NumColumns, matrix.EnumerateNonZeros());
        }

        /// <summary>
        /// Initializes a new instance of <see cref="DokSymmetric"/> with the specified matrix dimensions and the 
        /// non-zero entries defined by the provided symmetric pattern.
        /// </summary>
        /// <param name="numRows">The number of rows/columns of the symmetric matrix to build.</param>
        /// <param name="nonZeroEntries">The non-zero entries of the symmetric matrix to build.</param>
        public static DokSymmetric CreateFromSparseSymmetricPattern(int order,
            IEnumerable<(int row, int col, double value)> nonZeroEntries)
        {
            DokSymmetric dok = CreateEmpty(order);
            foreach (var (row, col, val) in nonZeroEntries) dok.columns[col].Add(row, val);
            return dok;
        }

        /// <summary>
        /// See <see cref="IMatrixBuilder.AddToEntry(int, int, double)"/>
        /// </summary>
        public void AddToEntry(int rowIdx, int colIdx, double value)
        {
            ProcessIndices(ref rowIdx, ref colIdx);
            if (columns[colIdx].TryGetValue(rowIdx, out double oldValue))
            {
                columns[colIdx][rowIdx] = value + oldValue;
            }
            else columns[colIdx][rowIdx] = value;
            //The Dictionary columns[rowIdx] is indexed twice in both cases. Is it possible to only index it once?
        }

        /// <summary>
        /// See <see cref="ISymmetricMatrixBuilder.AddSubmatrix(IIndexable2D, IReadOnlyDictionary{int, int}, 
        /// IReadOnlyDictionary{int, int})"/>.
        /// </summary>
        public void AddSubmatrix(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols)
        {
            foreach (var colPair in subColsToGlobalCols)
            {
                int subCol = colPair.Key;
                foreach (var rowPair in subRowsToGlobalRows)
                {
                    int subRow = rowPair.Key;
                    int globalRow, globalCol;
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
        /// See <see cref="ISymmetricMatrixBuilder.AddSubmatrixSymmetric(IIndexable2D, IReadOnlyDictionary{int, int})"/>.
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
                    if (globalRow <= globalCol)
                    {
                        double subValue = subMatrix[subRow, subCol];
                        if (columns[globalCol].TryGetValue(globalRow, out double oldGlobalValue))
                        {
                            columns[globalCol][globalRow] = subValue + oldGlobalValue;
                        }
                        else columns[globalCol][globalRow] = subValue;
                    }
                }
            }
        }

        /// <summary>
        /// See <see cref="ISymmetricMatrixBuilder.AddSubmatrixSymmetric(IIndexable2D, int[], int[])"/>.
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
                    if (globalRow <= globalCol)
                    {
                        double subValue = subMatrix[subRow, subCol];
                        if (columns[globalCol].TryGetValue(globalRow, out double oldGlobalValue))
                        {
                            columns[globalCol][globalRow] = subValue + oldGlobalValue;
                        }
                        else columns[globalCol][globalRow] = subValue;
                    }
                }
            }
        }

        /// <summary>
        /// See <see cref="ISymmetricMatrixBuilder.AddSubmatrixToLowerTriangle(IIndexable2D, IReadOnlyDictionary{int, int}, 
        /// IReadOnlyDictionary{int, int})"/>
        /// </summary>
        public void AddSubmatrixToLowerTriangle(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols)
        {
            foreach (var rowPair in subRowsToGlobalRows) // Transpose(Col major ordering) = Row major ordering
            {
                int subRow = rowPair.Key;
                int globalRow = rowPair.Value;
                foreach (var colPair in subColsToGlobalCols)
                {
                    int subCol = colPair.Key;
                    int globalCol = colPair.Value;
                    double subValue = subMatrix[subRow, subCol];
                    // Transpose the access on the global matrix only. 
                    // This is equivalent to calling the indexer for each lower triangle entry, but without the explicit swap.
                    if (columns[globalRow].TryGetValue(globalCol, out double oldGlobalValue))
                    {
                        columns[globalRow][globalCol] = subValue + oldGlobalValue;
                    }
                    else columns[globalRow][globalCol] = subValue;
                }
            }
        }

        /// <summary>
        /// See <see cref="ISymmetricMatrixBuilder.AddSubmatrixToUpperTriangle(IIndexable2D, IReadOnlyDictionary{int, int}, 
        /// IReadOnlyDictionary{int, int})"/>
        /// </summary>
        public void AddSubmatrixToUpperTriangle(IIndexable2D subMatrix,
            IReadOnlyDictionary<int, int> subRowsToGlobalRows, IReadOnlyDictionary<int, int> subColsToGlobalCols)
        {
            foreach (var colPair in subColsToGlobalCols) // Col major ordering
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

        public (double[] values, int[] diagOffsets) BuildSkylineArrays()
        {
            //TODO: It may be worth it to create a DOK specifically for skyline. As you add entries to it, it will cache the
            //      min row or the height for each column. Alternatively it could use SortedDictionary, but that is slow for 
            //      global matrix assembly

            var diagOffsets = new int[order + 1];
            int nnz = 0;
            for (int j = 0; j < order; ++j)
            {
                diagOffsets[j] = nnz;
                nnz += j - columns[j].Keys.Min() + 1; // Height of column j including the diagonal entry
            }
            if (nnz == 0) throw new EmptyMatrixBuilderException("Cannot build symmetric CSC arrays from a DOK with nnz = 0.");
            diagOffsets[order] = nnz; //The last Skyline entry has diagOffset = nnz.

            var values = new double[nnz];
            for (int j = 0; j < order; ++j)
            {
                int temp = diagOffsets[j] + j;
                foreach (var rowVal in columns[j])
                {
                    // offset = diagOffsets[j] + j - row
                    int offset = temp - rowVal.Key;
                    values[offset] = rowVal.Value;
                }
            }
            return (values, diagOffsets);
        }

        public SkylineMatrix BuildSkylineMatrix()
        {
            (double[] values, int[] diagOffsets) = BuildSkylineArrays();
            return SkylineMatrix.CreateFromArrays(NumColumns, values, diagOffsets, false, false);
        }

        /// <summary>
        /// Creates the values and indexing arrays in CSC storage format of the upper triangle of the current symmetric matrix. 
        /// This method should be called after fully defining the matrix in <see cref="DokSymmetric"/> format.
        /// </summary>
        /// <param name="sortColsOfEachCol">True to sort the column indices of the CSC matrix between colOffsets[j] and 
        ///     colOffsets[j+1] in ascending order. False to leave them unordered. Ordered rows might result in better  
        ///     performance during multiplications or they might be required by 3rd party libraries. Conversely, leaving them 
        ///     unordered will be faster during creation of the CSC matrix.</param>
        /// <exception cref="EmptyMatrixBuilderException">Thrown if no non-zero entries have been defined yet.</exception>
        public (double[] values, int[] rowIndices, int[] columnOffsets) BuildSymmetricCscArrays(bool sortRowsOfEachCol)
        {
            int[] colOffsets = new int[order + 1];
            int nnz = 0;
            for (int j = 0; j < order; ++j)
            {
                colOffsets[j] = nnz;
                nnz += columns[j].Count;
            }
            if (nnz == 0) throw new EmptyMatrixBuilderException("Cannot build symmetric CSC arrays from a DOK with nnz = 0.");
            colOffsets[order] = nnz; //The last CSC entry has colOffset = nnz.

            int[] rowIndices = new int[nnz];
            double[] values = new double[nnz];
            int counter = 0;
            if (sortRowsOfEachCol)
            {
                for (int j = 0; j < order; ++j)
                { //TODO: perhaps passing the dictionarirs to new SortedDictionary() is faster than LINQ
                    foreach (var rowVal in columns[j].OrderBy(pair => pair.Key)) 
                    {
                        rowIndices[counter] = rowVal.Key;
                        values[counter] = rowVal.Value;
                        ++counter;
                    }
                }
            }
            else
            {
                for (int j = 0; j < order; ++j)
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
        /// Initializes a <see cref="SymmetricCscMatrix"/> representation of the current matrix. This method should be 
        /// called after fully defining the matrix in <see cref="DokSymmetric"/> format.
        /// </summary>
        /// <param name="sortColsOfEachCol">True to sort the column indices of the CSC matrix between colOffsets[j] and 
        ///     colOffsets[j+1] in ascending order. False to leave them unordered. Ordered rows might result in better  
        ///     performance during multiplications or they might be required by 3rd party libraries. Conversely, leaving them 
        ///     unordered will be faster during creation of the CSC matrix.</param>
        /// <exception cref="EmptyMatrixBuilderException">Thrown if no non-zero entries have been defined yet.</exception>
        public SymmetricCscMatrix BuildSymmetricCscMatrix(bool sortRowsOfEachCol)
        {
            (double[] values, int[] rowIndices, int[] colOffsets) = BuildSymmetricCscArrays(sortRowsOfEachCol);
            return SymmetricCscMatrix.CreateFromArrays(NumColumns, values, rowIndices, colOffsets, false);
        }

        /// <summary>
        /// Frees all memory held by this <see cref="DokSymmetric"/> instance. Afterwards this object cannot be reused. 
        /// Therefore this method should be called to save up space after building the matrix or its internal arrays.
        /// </summary>
        public void Clear() => columns = null;

        /// <summary>
        /// See <see cref="ISparseMatrix.CountNonZeros"/>.
        /// </summary>
        public int CountNonZeros()
        {
            int count = 0;
            for (int j = 0; j < order; ++j)
            {
                //count += (columns[j].Count - 1) * 2 + 1; // Faster, but doesn't work if the diagonal is not present.
                foreach (var rowVal in columns[j])
                {
                    if (rowVal.Key == j) ++count;
                    else count += 2; //Each upper triangle entries has a corresponding lower triangle entry.
                }
            }
            return count;
        }

        /// <summary>
        /// See <see cref="ISparseSymmetricMatrix.CountNonZerosUpper"/>.
        /// </summary>
        public int CountNonZerosUpper()
        {
            int count = 0;
            for (int j = 0; j < order; ++j) count += columns[j].Count;
            return count;
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.EnumerateNonZeros"/>.
        /// </summary>
        public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
        {
            for (int j = 0; j < order; ++j)
            {
                foreach (var rowVal in columns[j])
                {
                    if (rowVal.Key == j) yield return (rowVal.Key, j, rowVal.Value);
                    else //Each upper triangle entries has a corresponding lower triangle entry.
                    {
                        yield return (rowVal.Key, j, rowVal.Value);
                        yield return (j, rowVal.Key, rowVal.Value);
                    }
                }
            }
        }

        /// <summary>
        /// See <see cref="ISparseSymmetricMatrix.EnumerateNonZerosUpper"/>.
        /// </summary>
        public IEnumerable<(int row, int col, double value)> EnumerateNonZerosUpper()
        {
            for (int j = 0; j < order; ++j)
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
        /// Returns the column with index = <paramref name="colIdx"/> as a sparse vector. Note that the length of the returned 
        /// vector is equal to this.<see cref="NumRows"/>. Since the matrix is symmetric, this method also works for getting
        /// the row with index = <paramref name="colIdx"/>.
        /// </summary>
        /// <param name="colIdx">The index of the column to return. Constraints: Column <paramref name="colIdx"/> must be stored
        ///     and 0 &lt;= <paramref name="colIdx"/> &lt; this.<see cref="NumColumns"/>.</param>
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="colIdx"/> or <paramref name="tabooRows"/> 
        ///     violate the described constraints.</exception>
        public SparseVector GetColumn(int colIdx)
        { //TODO: implement another data structure that only holds some sparse columns, but is efficient in returning them whole
            // The super-diagonal part is readily available.
            var wholeColumn = new SortedDictionary<int, double>(columns[colIdx]);

            // The sub-diagonal part of the column is stored as the super-diagonal part of the row with the same index
            // It must be built by searching all subsequent columns, which is inefficient
            for (int j = colIdx + 1; j < order; ++j)
            {
                bool isNonZero = columns[j].TryGetValue(colIdx, out double value);
                if (isNonZero) wholeColumn.Add(j, value);
            }
            return SparseVector.CreateFromDictionary(order, wholeColumn);
        }

        /// <summary>
        /// Returns the column with index = <paramref name="colIdx"/> as a vector. However, the entries with row index that 
        /// belongs in <paramref name="tabooRows"/> will be set to 0. More accurately, they will not be included in the 
        /// sparsity pattern of the returned <see cref="SparseVector"/>. Note that the length of the returned vector is 
        /// equal to this.<see cref="NumRows"/>. Since the matrix is symmetric, this method also works for getting the row
        /// with index = <paramref name="colIdx"/>.
        /// </summary>
        /// <param name="colIdx">The index of the column to return. Constraints: Column <paramref name="colIdx"/> must be stored
        ///     and 0 &lt;= <paramref name="colIdx"/> &lt; this.<see cref="NumColumns"/>.</param>
        /// <param name="tabooRows">The entries of the returned column vector at the indices specified by 
        ///     <paramref name="tabooRows"/> will be equal to 0. Constraints: foreach rowIdx in <paramref name="tabooRows"/>:
        ///     0 &lt;= rowIdx &lt; this.<see cref="NumRows"/>.</param>
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="colIdx"/> or <paramref name="tabooRows"/> 
        ///     violate the described constraints.</exception>
        public SparseVector GetColumnWithoutRows(int colIdx, ISet<int> tabooRows)
        { //TODO: it's not really slicing, as the rejected rows are set to 0
            // The super-diagonal part is straightforward
            var wholeColumn = new SortedDictionary<int, double>();
            foreach (var rowVal in columns[colIdx])
            {
                if (!tabooRows.Contains(rowVal.Key)) wholeColumn.Add(rowVal.Key, rowVal.Value);
            }

            // The sub-diagonal part of the column is stored as the super-diagonal part of the row with the same index
            // It must be accessed by searching all subsequent columns, which is inefficient
            for (int i = colIdx + 1; i < order; ++i)
            {
                if (!tabooRows.Contains(i))
                {
                    bool isNonZero = columns[i].TryGetValue(colIdx, out double value);
                    if (isNonZero) wholeColumn.Add(i, value);
                }
            }
            return SparseVector.CreateFromDictionary(order, wholeColumn);
        }

        /// <summary>
        /// Returns a <see cref="Vector"/> with the entries of the matrix's main diagonal and the index of the first zero entry.
        /// If there are no zero entries, -1 is returned as the index.
        /// </summary>
        public (Vector diagonal, int firstZeroIdx) GetDiagonal()
        {
            (double[] diagonal, int firstZeroIdx) = GetDiagonalAsArray();
            return (Vector.CreateFromArray(diagonal, false), firstZeroIdx);
        }

        /// <summary>
        /// Returns an array with the entries of the matrix's main diagonal and the index of the first zero entry.
        /// If there are no zero entries, -1 is returned as the index.
        /// </summary>
        public (double[] diagonal, int firstZeroIdx) GetDiagonalAsArray()
        {
            double[] diag = new double[order];
            int firstZeroIdx = -1;
            for (int j = 0; j < order; ++j)
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
        { //Perhaps there should be a dedicated symmetric CSC format, identical to CSC.
            (double[] values, int[] rowIndices, int[] colOffsets) = BuildSymmetricCscArrays(false);
            var format = new SparseFormat();
            format.RawValuesTitle = "Values";
            format.RawValuesArray = values;
            format.RawIndexArrays.Add("Row indices", rowIndices);
            format.RawIndexArrays.Add("Column offsets", colOffsets);
            return format;
        }

        /// <summary>
        /// Creates a new <see cref="DokSymmetric"/> by copying only the entries whose row and column index belongs to 
        /// <paramref name="rowsColsToKeep"/>.
        /// </summary>
        /// <param name="rowsColsToKeep">
        /// Only the entries (i, j) where i and j belongs to <paramref name="rowsColsToKeep"/> will be copied to the new matrix.
        /// The order of the indices is important. E.g. if <paramref name="rowsColsToKeep"/> = [1, 0], then the resulting 
        /// submatrix will be {{ this[1,1], this[1,0] }, { this[0,1], this[0,0] }}
        /// </param>
        public DokSymmetric GetSubmatrix(int[] rowsColsToKeep)
        {
            //TODO: work with the stored dictionaries to speed this up. However keep in mind that if keep = {2, 1} then some 
            //      transposition is necessary. 
            int subOrder = rowsColsToKeep.Length;
            var submatrix = DokSymmetric.CreateEmpty(subOrder);
            for (int subJ = 0; subJ < subOrder; ++subJ)
            {
                int thisJ = rowsColsToKeep[subJ];
                for (int subI = 0; subI <= subJ; ++subI)
                {
                    int thisI = rowsColsToKeep[subI];
                    double value = this[thisI, thisJ]; //TODO: This should be done explicitly
                    if (value != 0.0) submatrix[subI, subJ] = value;
                }
            }
            return submatrix;
        }

        /// <summary>
        /// Sets the column with index <paramref name="colIdx"/> to be equal to the provided <paramref name="newColumn"/>.
        /// Since the matrix is symmetric, this method also works for modifying the row with index = <paramref name="colIdx"/>.
        /// </summary>
        /// <param name="colIdx">The index of the column to modify. Constraints:
        ///     and 0 &lt;= <paramref name="colIdx"/> &lt; this.<see cref="NumColumns"/>.</param>
        /// <param name="newColumn">The new values that column <paramref name="colIdx"/> will be set to. Constraints:
        ///     <paramref name="newColumn"/>.<see cref="IIndexable1D.Length"/> == this.<see cref="NumRows"/>.</param>
        public void SetColumn(int colIdx, SparseVector newColumn)
        {
            SetColumnToZero(colIdx); // First remove everything

            int[] rowIndices = newColumn.RawIndices;
            double[] values = newColumn.RawValues;

            // The super-diagonal part is straightforward
            int t = 0;
            for (; ; ++t)
            {
                int rowIdx = rowIndices[t];
                if (rowIdx > colIdx) break;
                columns[colIdx].Add(rowIdx, values[t]);
            }

            // The sub-diagonal part of the column is stored as the super-diagonal part of the row with the same index
            // It must be accessed by searching the relevant subsequent columns
            for (; t < rowIndices.Length; ++t)
            {
                int rowIdx = rowIndices[t];
                columns[rowIdx].Add(colIdx, values[t]); // It should have been emptied when it was set to identity
            }
        }

        /// <summary>
        /// Sets the column with index <paramref name="colIdx"/> to have 1.0 at the main diagonal entry and 0.0 everywhere else.
        /// Since the matrix is symmetric, this method also works for modifying the row with index = <paramref name="colIdx"/>.
        /// </summary>
        /// <param name="colIdx">The index of the column to modify. Constraints:
        ///     and 0 &lt;= <paramref name="colIdx"/> &lt; this.<see cref="NumColumns"/>.</param>
        public void SetColumnToIdentity(int colIdx)
        {
            SetColumnToZero(colIdx); // First remove everything
            columns[colIdx].Add(colIdx, 1.0);
        }

        /// <summary>
        /// Sets the column with index <paramref name="colIdx"/> to have 0.0 at all entries.
        /// Since the matrix is symmetric, this method also works for modifying the row with index = <paramref name="colIdx"/>.
        /// </summary>
        /// <param name="colIdx">The index of the column to modify. Constraints:
        ///     and 0 &lt;= <paramref name="colIdx"/> &lt; this.<see cref="NumColumns"/>.</param>
        public void SetColumnToZero(int colIdx)
        {
            // The super-diagonal part of the column is straightforward
            columns[colIdx].Clear();

            // The sub-diagonal part of the column is stored as the super-diagonal part of the row with the same index
            // It must be accessed by searching all subsequent columns, which is inefficient.
            for (int j = colIdx + 1; j < order; ++j)
            {
                columns[j].Remove(colIdx); // If it is 0, then nothing will happen.
            }
        }

        /// <summary>
        /// Like <see cref="this[int, int]"/>, but will not check or transpose the entry.
        /// </summary>
        /// <param name="rowIdx">The row index: 0 &lt;= <paramref name="rowIdx"/> &lt; <see cref="NumRows"/>.</param>
        /// <param name="colIdx">The column index: 0 &lt;= <paramref name="colIdx"/> &lt; <see cref="NumColumns"/>.</param>
        /// <param name="value">The value to set.</param>
        public void SetEntryUpper(int rowIdx, int colIdx, double value) => columns[colIdx][rowIdx] = value;

        /// <summary>
        /// Sets all entries that are on the main diagonal and have not been explicitly modified by a previous method or 
        /// constructor, and thus are equal to 0, to 1.
        /// </summary>
        public void SetStructuralZeroDiagonalEntriesToUnity()
        {
            for (int d = 0; d < NumColumns; ++d)
            {
                bool structuralZero = !columns[d].ContainsKey(d);
                if (structuralZero) columns[d][d] = 1.0;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void ProcessIndices(ref int rowIdx, ref int colIdx)
        { //Perhaps this should be manually inlined. Testing needed.
            if (rowIdx > colIdx)
            {
                int swap = rowIdx;
                rowIdx = colIdx;
                colIdx = swap;
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
        private void AddSubmatrixSymmetricOLD(IIndexable2D elementMatrix, int[] elementDofs, int[] globalDofs) //TODO: this should be reworked
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
                for (int i = 0; i <= j; ++i)
                {
                    int elementRow = elementDofs[j];
                    int globalRow = globalDofs[j];
                    double newGlobalValue = elementMatrix[elementRow, elementCol];
                    if (columns[globalCol].TryGetValue(globalRow, out double oldGlobalValue))
                    {
                        newGlobalValue += oldGlobalValue;
                    }
                    columns[globalCol][globalRow] = newGlobalValue;
                }
            }
        }
    }
}
