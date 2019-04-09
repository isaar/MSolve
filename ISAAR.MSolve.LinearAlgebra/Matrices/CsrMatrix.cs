using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using static ISAAR.MSolve.LinearAlgebra.LibrarySettings;

//TODO: try to make general versions of row major and col major multiplication. The lhs matrix/vector will be supplied by the 
//      caller, depending on if it is a matrix, a transposed matrix, a vector, etc. Compare its performance with the verbose 
//      code and also try inlining. C preprocessor macros would be actually useful here. Otherwise, move all that boilerplate 
//      code to a CSRStrategies static class.
//TODO: In matrix-matrix/vector multiplications: perhaps I should work with a column major array directly instead of an output   
//      Matrix and an array instead of an output Vector.
//TODO: perhaps optimizations if (other is Matrix) are needed, to directly index into its raw col major array.
//      The access paterns are always the same
//TODO: have a field: bool AreRowsSorted. Then it can be inspected to decide on more efficient algorithms (e.g. binary search 
//      of a row to find a column index.
namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Sparse matrix stored in Compressed Sparse Rows format (3-array version). The CSR format is optimized for matrix-vector 
    /// and matrix-matrix multiplications, where the CSR matrix is on the left untransposed or on the right transposed. The other
    /// multiplicationss are more efficient using <see cref="CscMatrix"/>. To build a <see cref="CsrMatrix"/> conveniently, 
    /// use <see cref="Builders.DokRowMajor"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CsrMatrix : IMatrix, ISparseMatrix
    {
        private const int zeroEntryOffset = -1;

        private readonly double[] values;
        private readonly int[] colIndices;
        private readonly int[] rowOffsets;

        private CsrMatrix(int numRows, int numCols, double[] values, int[] colIndices, int[] rowOffsets)
        {
            this.values = values;
            this.colIndices = colIndices;
            this.rowOffsets = rowOffsets;
            this.NumRows = numRows;
            this.NumColumns = numCols;
        }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of non zero entries of the matrix.
        /// </summary>
        public int NumNonZeros => rowOffsets[rowOffsets.Length - 1];

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get; }

        /// <summary>
        /// The internal array that stores the non-zero entries of the matrix. The non-zero entries of each row are consecutive.
        /// Its length is equal to the number of non-zero entries. 
        /// It should only be used for passing the raw array to linear algebra libraries.
        /// </summary>
        internal double[] RawValues => values;

        /// <summary>
        /// The internal array that stores the column indices of the non-zero entries in <see cref="RawValues"/>.
        /// Its length is equal to the number of non-zero entries. 
        /// It should only be used for passing the raw array to linear algebra libraries.
        /// </summary>
        internal int[] RawColIndices => colIndices;

        /// <summary>
        /// The internal array that stores the index into the arrays <see cref="RawValues"/> and <see cref="RawColIndices"/> of  
        /// the first entry of each row. Its length is equal to <paramref name="NumRows"/> + 1. The last entry is the number of 
        /// non-zero entries, which must be equal to <see cref="RawValues"/>.Length == <see cref="RawColIndices"/>.Length.
        /// It should only be used for passing the raw array to linear algebra libraries.
        /// </summary>
        internal int[] RawRowOffsets => rowOffsets;

        /// <summary>
        /// See <see cref="IIndexable2D.this[int, int]"/>.
        /// </summary>
        public double this[int rowIdx, int colIdx]
        {
            get
            {
                int entryOffset = FindOffsetOf(rowIdx, colIdx);
                if (entryOffset == zeroEntryOffset) return 0.0;
                else return values[entryOffset];
            }
        }

        /// <summary>
        /// Initializes a new <see cref="CsrMatrix"/> with the specified dimensions and the provided arrays 
        /// (<paramref name="values"/>, <paramref name="colIndices"/> and <paramref name="rowOffsets"/>) as its internal data.
        /// </summary>
        /// <param name="numRows">The number of rows of the new matrix.</param>
        /// <param name="numCols">The number of columns of the new matrix.</param>
        /// <param name="values">Array that contains the non-zero entries. It must have the same length 
        ///     as <paramref name="colIndices"/>. The non-zero entries of each row must appear consecutively in 
        ///     <paramref name="values"/>. They can also be sorted in increasing order of their column indices, which speeds up
        ///     subsequent operations.</param>
        /// <param name="colIndices">Array that contains the column indices of the non-zero entries. It must have the same 
        ///     length as <paramref name="values"/>. There is an 1 to 1 matching between these two arrays: 
        ///     <paramref name="colIndices"/>[i] is the column index of the entry <paramref name="values"/>[i]. Also:
        ///     0 &lt;= <paramref name="colIndices"/>[i] &lt; <paramref name="numCols"/>.</param>
        /// <param name="rowOffsets">Array that contains the index of the first entry of each row into the arrays 
        ///     <paramref name="values"/> and <paramref name="colIndices"/>. Its length is <paramref name="numRows"/> + 1. The 
        ///     last entry is the number of non-zero entries, which must be equal to the length of <paramref name="values"/> 
        ///     and <paramref name="colIndices"/>.</param>
        /// <param name="checkInput">If true, the provided arrays will be checked to make sure they are valid CSR arrays, which 
        ///     is safer. If false, no such check will take place, which is faster.</param>
        public static CsrMatrix CreateFromArrays(int numRows, int numCols, double[] values, int[] colIndices, int[] rowOffsets, 
            bool checkInput)
        {
            if (checkInput)
            {
                if (rowOffsets.Length != numRows + 1)
                {
                    throw new ArgumentException("The length of the CSR row offsets array must be equal to the number of rows + 1"
                        + ", but was " + rowOffsets.Length);
                }
                if (values.Length != colIndices.Length)
                {
                    throw new ArgumentException("The length of the CSR values and column indices arrays must be equal (and equal"
                        + $" to the number of non zero entries), but were {values.Length} and {colIndices.Length} respectively");
                }
                if (rowOffsets[0] != 0)
                {
                    throw new ArgumentException("The first entry of the CSR column offsets array must be 0, but was "
                        + rowOffsets[0]);
                }
                if (rowOffsets[rowOffsets.Length-1] != values.Length)
                {
                    throw new ArgumentException("The last entry of the CSR row offsets array must be equal to the number of non"
                        + " zero entries, but was " + rowOffsets[rowOffsets.Length - 1]);
                }
            }
            return new CsrMatrix(numRows, numCols, values, colIndices, rowOffsets);
        }

        #region operators (use extension operators when they become available)
        /// <summary>
        /// Performs the matrix-vector multiplication: result = <paramref name="matrixLeft"/> * <paramref name="vectorRight"/>.
        /// If <paramref name="matrixLeft"/> is m1-by-n1 and <paramref name="vectorRight"/> has length = n2, then n1 must be 
        /// equal to n2. The result will be a vector with length = m1, written to a new <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="matrixLeft">The <see cref="CsrMatrix"/> operand on the left.</param>
        /// <param name="vectorRight">The <see cref="Vector"/> operand on the right. It can be considered as a column 
        ///     vector.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrixLeft"/>.<see cref="NumColumns"/> is 
        ///     different than <paramref name="vectorRight"/>.<see cref="Vector.Length"/>.</exception>
        public static Vector operator *(CsrMatrix matrixLeft, Vector vectorRight)
            => matrixLeft.Multiply(vectorRight, false);
        #endregion

        /// <summary>
        /// See <see cref="IMatrixView.Axpy(IMatrixView, double)"/>.
        /// </summary>
        public IMatrix Axpy(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CsrMatrix otherCSR) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherCSR))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    Array.Copy(this.values, resultValues, values.Length);
                    Blas.Daxpy(values.Length, otherCoefficient, otherCSR.values, 0, 1, resultValues, 0, 1);
                    return new CsrMatrix(NumRows, NumColumns, resultValues, this.colIndices, this.rowOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.LinearCombination(this, 1.0, otherMatrix, otherCoefficient);
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// result[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
        /// The resulting matrix is written to a new <see cref="CsrMatrix"/> and then returned.
        /// </summary>
        /// <param name="otherMatrix">A matrix with the same <see cref="NumRows"/> and <see cref="NumColumns"/> as this 
        ///     <see cref="CsrMatrix"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     <see cref="NumRows"/> or <see cref="NumColumns"/> than this instance.</exception>
        public CsrMatrix Axpy(CsrMatrix otherMatrix, double otherCoefficient)
        {
            // Conceptually it is not wrong to so this, even if the indexers are different, but how would I implement it.
            if (!HasSameIndexer(otherMatrix)) 
            {
                throw new SparsityPatternModifiedException("Only allowed if the indexing arrays are the same");
            }
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] resultValues = new double[values.Length];
            Array.Copy(this.values, resultValues, values.Length);
            Blas.Daxpy(values.Length, otherCoefficient, otherMatrix.values, 0, 1, resultValues, 0, 1);
            // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
            return new CsrMatrix(NumRows, NumColumns, resultValues, this.colIndices, this.rowOffsets);
        }

        /// <summary>
        /// See <see cref="IMatrix.AxpyIntoThis(IMatrixView, double)"/>.
        /// </summary>
        public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CsrMatrix casted) AxpyIntoThis(casted, otherCoefficient);
            else throw new SparsityPatternModifiedException(
                 "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// this[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
        /// The resulting matrix overwrites the entries of this <see cref="CsrMatrix"/> instance.
        /// </summary>
        /// <param name="otherMatrix">A matrix with the same indexing arrays as this <see cref="CsrMatrix"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="SparsityPatternModifiedException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     indexing arrays than this instance.</exception>
        public void AxpyIntoThis(CsrMatrix otherMatrix, double otherCoefficient)
        {
            //Preconditions.CheckSameMatrixDimensions(this, other); // no need if the indexing arrays are the same
            if (!HasSameIndexer(otherMatrix))
            {
                throw new SparsityPatternModifiedException("Only allowed if the indexing arrays are the same");
            }
            Blas.Daxpy(values.Length, otherCoefficient, otherMatrix.values, 0, 1, this.values, 0, 1);
        }

        /// <summary>
        /// See <see cref="IMatrixView.CopyToFullMatrix()"/>
        /// </summary>
        public Matrix CopyToFullMatrix()
        {
            Matrix fullMatrix = Matrix.CreateZero(this.NumRows, this.NumColumns);
            for (int i = 0; i < NumRows; ++i) //Row major order
            {
                int rowStart = rowOffsets[i]; //inclusive
                int rowEnd = rowOffsets[i + 1]; //exclusive
                for (int k = rowStart; k < rowEnd; ++k)
                {
                    fullMatrix[i, colIndices[k]] = values[k];
                }
            }
            return fullMatrix;
        }

        /// <summary>
        /// See <see cref="IMatrix.Clear"/>.
        /// </summary>
        public void Clear() => Array.Clear(values, 0, values.Length);

        /// <summary>
        /// See <see cref="IMatrixView.Copy(bool)"/>.
        /// </summary>
        IMatrix IMatrixView.Copy(bool copyIndexingData) => Copy(copyIndexingData);

        /// <summary>
        /// Copies the entries of this matrix.
        /// </summary>
        /// <param name="copyIndexingData">
        /// If true, all data of this object will be copied. If false, only the array containing the values of the stored 
        /// matrix entries will be copied. The new matrix will reference the same indexing arrays as this one.
        /// </param>
        public CsrMatrix Copy(bool copyIndexingData)
        {
            var valuesCopy = new double[this.values.Length];
            Array.Copy(this.values, valuesCopy, this.values.Length);

            if (!copyIndexingData) return new CsrMatrix(NumRows, NumColumns, valuesCopy, this.colIndices, this.rowOffsets);
            else
            {
                var colIndicesCopy = new int[this.colIndices.Length];
                Array.Copy(this.colIndices, colIndicesCopy, this.colIndices.Length);
                var rowOffsetsCopy = new int[this.rowOffsets.Length];
                Array.Copy(this.rowOffsets, rowOffsetsCopy, this.rowOffsets.Length);
                return new CsrMatrix(NumRows, NumColumns, valuesCopy, colIndicesCopy, rowOffsetsCopy);
            }
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.CountNonZeros"/>
        /// </summary>
        public int CountNonZeros() => values.Length;

        /// <summary>
        /// See <see cref="IMatrixView.DoEntrywise(IMatrixView, Func{double, double, double})"/>.
        /// </summary>
        public IMatrix DoEntrywise(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is CsrMatrix otherCSR) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherCSR))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    for (int i = 0; i < values.Length; ++i)
                    {
                        resultValues[i] = binaryOperation(this.values[i], otherCSR.values[i]);
                    }
                    return new CsrMatrix(NumRows, NumColumns, resultValues, colIndices, rowOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.DoEntrywise(this, other, binaryOperation);
        }

        /// <summary>
        /// See <see cref="IMatrix.DoEntrywiseIntoThis(IMatrixView, Func{double, double, double})"/>.
        /// </summary>
        public void DoEntrywiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is CsrMatrix casted) DoEntrywiseIntoThis(casted, binaryOperation);
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// this[i, j] = <paramref name="binaryOperation"/>(this[i,j], <paramref name="matrix"/>[i, j]). 
        /// The resulting matrix overwrites the entries of this <see cref="CsrMatrix"/> instance.
        /// </summary>
        /// <param name="matrix">A matrix with the same indexing arrays as this <see cref="CsrMatrix"/> instance.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        /// <exception cref="SparsityPatternModifiedException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     indexing arrays than this instance.</exception>
        public void DoEntrywiseIntoThis(CsrMatrix other, Func<double, double, double> binaryOperation)
        {
            //Preconditions.CheckSameMatrixDimensions(this, other); // no need if the indexing arrays are the same
            if (!HasSameIndexer(other))
            {
                throw new SparsityPatternModifiedException("Only allowed if the indexing arrays are the same");
            }
            for (int i = 0; i < values.Length; ++i) this.values[i] = binaryOperation(this.values[i], other.values[i]);
        }

        /// <summary>
        /// See <see cref="IMatrixView.DoToAllEntries(Func{double, double})"/>.
        /// </summary>
        IMatrix IMatrixView.DoToAllEntries(Func<double, double> unaryOperation)
        {
            // Only apply the operation on non zero entries
            double[] newValues = new double[values.Length];
            for (int i = 0; i < values.Length; ++i) newValues[i] = unaryOperation(values[i]);

            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0)) // The same sparsity pattern can be used.
            {
                // Copy the index arrays. TODO: See if we can use the same index arrays (e.g. if this class does not change them (it shouldn't))
                int[] colIndicesCopy = new int[colIndices.Length];
                Array.Copy(colIndices, colIndicesCopy, colIndices.Length);
                int[] rowOffsetsCopy = new int[rowOffsets.Length];
                Array.Copy(rowOffsets, rowOffsetsCopy, rowOffsets.Length);
                return new CsrMatrix(NumRows, NumColumns, newValues, colIndicesCopy, rowOffsetsCopy);
            }
            else // The sparsity is destroyed. Revert to a full matrix.
            {
                return new CsrMatrix(NumRows, NumColumns, newValues, colIndices, rowOffsets).CopyToFullMatrix();
            }
        }

        /// <summary>
        /// See <see cref="IMatrix.DoToAllEntriesIntoThis(Func{double, double})"/>.
        /// </summary>
        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0))
            {
                for (int i = 0; i < values.Length; ++i) values[i] = unaryOperation(values[i]);
            }
            else
            {
                throw new SparsityPatternModifiedException("This operation will change the sparsity pattern");
            }
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.EnumerateNonZeros"/>.
        /// </summary>
        public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
        {
            for (int i = 0; i < NumRows; ++i)
            {
                int rowStart = rowOffsets[i]; //inclusive
                int rowEnd = rowOffsets[i + 1]; //exclusive
                for (int k = rowStart; k < rowEnd; ++k)
                {
                    yield return (i, colIndices[k], values[k]);
                }
            }
        }

        /// <summary>
        /// See <see cref="IIndexable2D.Equals(IIndexable2D, double)"/>.
        /// </summary>
        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
        {
            if ((this.NumRows != other.NumRows) || (this.NumColumns != other.NumColumns)) return false;
            var comparer = new ValueComparer(1e-13);
            for (int i = 0; i < NumRows; ++i)
            {
                int rowStart = rowOffsets[i]; //inclusive
                int rowEnd = rowOffsets[i + 1]; //exclusive
                int previousCol = 0;
                for (int k = rowStart; k < rowEnd; ++k)
                {
                    int col = colIndices[k];
                    for (int j = previousCol; j < col; ++j) // zero entries between the stored ones
                    {
                        if (!comparer.AreEqual(0.0, other[i, j])) return false;
                    }
                    if (!comparer.AreEqual(values[k], other[i, col])) return false; // Non zero entry
                    previousCol = col + 1;
                }
            }
            return true; //At this point all entries have been checked and are equal
        }

        /// <summary>
        /// See <see cref="ISliceable2D.GetColumn(int)"/>.
        /// </summary>
        public Vector GetColumn(int colIndex)
        {
            Preconditions.CheckIndexCol(this, colIndex);
            double[] colVector = new double[NumRows];
            for (int i = 0; i < NumRows; ++i)
            {
                int entryOffset = FindOffsetOf(i, colIndex);
                if (entryOffset != zeroEntryOffset) colVector[i] = values[entryOffset];
            }
            return Vector.CreateFromArray(colVector, false);
        }

        /// <summary>
        /// See <see cref="ISliceable2D.GetRow(int)"/>.
        /// </summary>
        public Vector GetRow(int rowIndex)
        {
            Preconditions.CheckIndexRow(this, rowIndex);
            double[] rowVector = new double[NumColumns];
            for (int k = rowOffsets[rowIndex]; k < rowOffsets[rowIndex + 1]; ++k) rowVector[colIndices[k]] = values[k];
            return Vector.CreateFromArray(rowVector, false);
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.GetSparseFormat"/>.
        /// </summary>
        public SparseFormat GetSparseFormat()
        {
            var format = new SparseFormat();
            format.RawValuesTitle = "Values";
            format.RawValuesArray = values;
            format.RawIndexArrays.Add("Column indices", colIndices);
            format.RawIndexArrays.Add("Row offsets", rowOffsets);
            return format;
        }

        /// <summary>
        /// See <see cref="ISliceable2D.GetSubmatrix(int[], int[])"/>.
        /// </summary>
        public Matrix GetSubmatrix(int[] rowIndices, int[] colIndices)
            => DenseStrategies.GetSubmatrix(this, rowIndices, colIndices);

        /// <summary>
        /// See <see cref="ISliceable2D.GetSubmatrix(int, int, int, int)"/>.
        /// </summary>
        public Matrix GetSubmatrix(int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive)
            => DenseStrategies.GetSubmatrix(this, rowStartInclusive, rowEndExclusive, colStartInclusive, colEndExclusive);

        /// <summary>
        /// See <see cref="IMatrixView.LinearCombination(double, IMatrixView, double)"/>.
        /// </summary>
        public IMatrix LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CsrMatrix otherCSR) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherCSR))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    if (thisCoefficient == 1.0)
                    {
                        Array.Copy(this.values, resultValues, values.Length);
                        Blas.Daxpy(values.Length, otherCoefficient, otherCSR.values, 0, 1, this.values, 0, 1);
                    }
                    else if (otherCoefficient == 1.0)
                    {
                        Array.Copy(otherCSR.values, resultValues, values.Length);
                        Blas.Daxpy(values.Length, thisCoefficient, this.values, 0, 1, resultValues, 0, 1);
                    }
                    else
                    {
                        Array.Copy(this.values, resultValues, values.Length);
                        BlasExtensions.Daxpby(values.Length, otherCoefficient, otherCSR.values, 0, 1,
                            thisCoefficient, this.values, 0, 1);
                    }
                    return new CsrMatrix(NumRows, NumColumns, resultValues, this.colIndices, this.rowOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.LinearCombination(this, thisCoefficient, otherMatrix, otherCoefficient);
        }

        /// <summary>
        /// See <see cref="IMatrix.LinearCombinationIntoThis(double, IMatrixView, double)"/>.
        /// </summary>
        public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CsrMatrix casted) LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// this[i, j] = <paramref name="thisCoefficient"/> * this[i, j] 
        ///     + <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j]. 
        /// The resulting matrix overwrites the entries of this <see cref="CsrMatrix"/> instance.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this <see cref="Matrix"/>.</param>
        /// <param name="otherMatrix">A matrix with the same indexing arrays as this <see cref="CsrMatrix"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="SparsityPatternModifiedException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     indexing arrays than this instance.</exception>
        public void LinearCombinationIntoThis(double thisCoefficient, CsrMatrix otherMatrix, double otherCoefficient)
        {
            //Preconditions.CheckSameMatrixDimensions(this, other); // no need if the indexing arrays are the same
            if (!HasSameIndexer(otherMatrix))
            {
                throw new SparsityPatternModifiedException("Only allowed if the indexing arrays are the same");
            }
            if (thisCoefficient == 1.0)
            {
                Blas.Daxpy(values.Length, otherCoefficient, otherMatrix.values, 0, 1, this.values, 0, 1);
            }
            else
            {
                BlasExtensions.Daxpby(values.Length, otherCoefficient, otherMatrix.values, 0, 1, 
                    thisCoefficient, this.values, 0, 1);
            }
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyLeft(IMatrixView, bool, bool)"/>.
        /// </summary>
        public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            //TODO: To use BLAS for this too, we must accept row major matrices as output.
            if (transposeOther)
            {
                if (transposeThis)
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumRows, this.NumColumns);
                    var result = Matrix.CreateZero(other.NumColumns, this.NumRows);
                    CsrMultiplications.MatrixTransTimesCsrTrans(this.NumRows, values, rowOffsets, colIndices, other, result);
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumRows, this.NumRows);
                    var result = Matrix.CreateZero(other.NumColumns, this.NumColumns);
                    CsrMultiplications.MatrixTransTimesCsr(this.NumRows, values, rowOffsets, colIndices, other, result);
                    return result;
                }
            }
            else
            {
                if (transposeThis)
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumColumns, this.NumColumns);
                    var result = Matrix.CreateZero(other.NumRows, this.NumRows);
                    CsrMultiplications.MatrixTimesCsrTrans(this.NumRows, values, rowOffsets, colIndices, other, result);
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumColumns, this.NumRows);
                    var result = Matrix.CreateZero(other.NumRows, this.NumColumns);
                    CsrMultiplications.MatrixTimesCsr(this.NumRows, values, rowOffsets, colIndices, other, result);
                    return result;
                }
            }
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyRight(IMatrixView, bool, bool)"/>.
        /// </summary>
        public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            if (transposeOther)
            {
                if (transposeThis)
                {
                    Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumColumns);
                    var result = Matrix.CreateZero(this.NumColumns, other.NumRows);
                    CsrMultiplications.CsrTransTimesMatrixTrans(this.NumRows, values, rowOffsets, colIndices, other, result);
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumColumns);
                    var result = Matrix.CreateZero(this.NumRows, other.NumRows);
                    CsrMultiplications.CsrTimesMatrixTrans(this.NumRows, values, rowOffsets, colIndices, other, result);
                    return result;
                }
            }
            else
            {
                //TODO: perhaps I can use the left multiplications if the other matrix is also transposed
                if (other is Matrix dense) return MultiplyRight(dense, transposeThis);

                if (transposeThis)
                {
                    Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumRows);
                    var result = Matrix.CreateZero(this.NumColumns, other.NumColumns);
                    CsrMultiplications.CsrTransTimesMatrix(this.NumRows, values, rowOffsets, colIndices, other, result);
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumRows);
                    var result = Matrix.CreateZero(this.NumRows, other.NumColumns);
                    CsrMultiplications.CsrTimesMatrix(this.NumRows, values, rowOffsets, colIndices, other, result);
                    return result;
                }
            }
        }

        /// <summary>
        /// Performs the matrix-matrix multiplication: oper(this) * <paramref name="other"/>.
        /// </summary>
        /// <param name="other">
        /// A matrix such that the <see cref="IIndexable2D.NumRows"/> of <paramref name="other"/> are equal to the 
        /// <see cref="IIndexable2D.NumColumns"/> of oper(this).
        /// </param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="otherMatrix"/> has different <see cref="IIndexable2D.NumRows"/> than the 
        /// <see cref="IIndexable2D.NumColumns"/> of oper(this).
        /// </exception>
        public Matrix MultiplyRight(Matrix other, bool transposeThis)
        {
            int numRowsResult;
            if (transposeThis)
            {
                Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumRows);
                numRowsResult = this.NumColumns;
            }
            else
            {
                Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumRows);
                numRowsResult = this.NumRows;
            }

            var result = Matrix.CreateZero(numRowsResult, other.NumColumns);
            SparseBlas.Dcsrgemm(transposeThis, this.NumRows, other.NumColumns, this.NumColumns, values, rowOffsets, colIndices,
                other.RawData, result.RawData);
            return result;
        }

        /// <summary>
        /// See <see cref="IMatrixView.Multiply(IVectorView, bool)"/>.
        /// </summary>
        public IVector Multiply(IVectorView vector, bool transposeThis = false)
        {
            if (vector is Vector casted) return Multiply(casted, transposeThis);

            if (transposeThis)
            {
                var result = new double[NumColumns];
                Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
                CsrMultiplications.CsrTransTimesVector(NumRows, values, rowOffsets, colIndices, vector, result);
                return Vector.CreateFromArray(result, false);
            }
            else
            {
                var result = new double[NumRows];
                Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
                CsrMultiplications.CsrTimesVector(NumRows, values, rowOffsets, colIndices, vector, result);
                return Vector.CreateFromArray(result, false);
            }
        }

        /// <summary>
        /// Performs the matrix-vector multiplication: oper(this) * <paramref name="vector"/>.
        /// To multiply this * columnVector, set <paramref name="transposeThis"/> to false.
        /// To multiply rowVector * this, set <paramref name="transposeThis"/> to true.
        /// </summary>
        /// <param name="vector">A vector with <see cref="IIndexable1D.Length"/> being equal to the 
        ///     <see cref="IIndexable2D.NumColumns"/> of oper(this).</param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the <see cref="IIndexable1D.Length"/> of
        ///     <paramref name="vector"/> is different than the <see cref="NumColumns"/> of oper(this).</exception>
        public Vector Multiply(Vector vector, bool transposeThis = false)
        {
            //TODO: this performs redundant dimension checks, including checking the transposeThis flag.
            var result = Vector.CreateZero(transposeThis ? NumColumns : NumRows);
            MultiplyIntoResult(vector, result, transposeThis);
            return result;
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyIntoResult(IVectorView, IVector, bool)"/>.
        /// </summary>
        public void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis = false)
        {
            if ((lhsVector is Vector lhsDense) && (rhsVector is Vector rhsDense))
            {
                MultiplyIntoResult(lhsDense, rhsDense, transposeThis);
            }

            if (transposeThis)
            {
                Preconditions.CheckMultiplicationDimensions(NumRows, lhsVector.Length);
                Preconditions.CheckSystemSolutionDimensions(NumColumns, rhsVector.Length);
                CsrMultiplications.CsrTransTimesVector(NumRows, values, rowOffsets, colIndices, lhsVector, rhsVector);
            }
            else
            {
                Preconditions.CheckMultiplicationDimensions(NumColumns, lhsVector.Length);
                Preconditions.CheckSystemSolutionDimensions(NumRows, rhsVector.Length);
                CsrMultiplications.CsrTimesVector(NumRows, values, rowOffsets, colIndices, lhsVector, rhsVector);
            }
        }

        /// <summary>
        /// Performs the matrix-vector multiplication: <paramref name="rhsVector"/> = oper(this) * <paramref name="vector"/>.
        /// To multiply this * columnVector, set <paramref name="transposeThis"/> to false.
        /// To multiply rowVector * this, set <paramref name="transposeThis"/> to true.
        /// The resulting vector will overwrite the entries of <paramref name="rhsVector"/>.
        /// </summary>
        /// <param name="lhsVector">
        /// The vector that will be multiplied by this matrix. It sits on the left hand side of the equation y = oper(A) * x.
        /// Constraints: <paramref name="lhsVector"/>.<see cref="IIndexable1D.Length"/> 
        /// == oper(this).<see cref="IIndexable2D.NumColumns"/>.
        /// </param>
        /// <param name="rhsVector">
        /// The vector that will be overwritten by the result of the multiplication. It sits on the right hand side of the 
        /// equation y = oper(A) * x. Constraints: <paramref name="lhsVector"/>.<see cref="IIndexable1D.Length"/> 
        /// == oper(this).<see cref="IIndexable2D.NumRows"/>.
        /// </param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if the <see cref="IIndexable1D.Length"/> of <paramref name="lhsVector"/> or <paramref name="rhsVector"/> 
        /// violate the described contraints.
        /// </exception>
        public void MultiplyIntoResult(Vector lhsVector, Vector rhsVector, bool transposeThis = false)
        {
            if (transposeThis)
            {
                Preconditions.CheckMultiplicationDimensions(NumRows, lhsVector.Length);
                Preconditions.CheckSystemSolutionDimensions(NumColumns, rhsVector.Length);
            }
            else
            {
                Preconditions.CheckMultiplicationDimensions(NumColumns, lhsVector.Length);
                Preconditions.CheckSystemSolutionDimensions(NumRows, rhsVector.Length);
            }
            SparseBlas.Dcsrgemv(transposeThis, NumRows, NumColumns, values, rowOffsets, colIndices, 
                lhsVector.RawData, 0, rhsVector.RawData, 0);
        }

        /// <summary>
        /// Performs the matrix-subvector multiplication (Matlab notation): 
        /// <paramref name="result"/>[<paramref name="resultStart"/>, :] = 
        ///     this * <paramref name="vectorRight"/>[<paramref name="vectorStart"/>, :]. 
        /// The resulting vector overwrites the entries of <paramref name="result"/> starting from the entry with index 
        /// <paramref name="resultStart"/>.
        /// </summary>
        /// <param name="vectorRight">The vector that will be multiplied with this matrix. <paramref name="vectorRight"/> is on 
        ///     the right of the multiplication.</param>
        /// <param name="vectorStart">The index of the first entry of <paramref name="vectorRight"/> that will be 
        ///     multiplied. Constraints: <paramref name="vectorStart"/> + this.<see cref="NumColumns"/> &lt;= 
        ///     <paramref name="vectorRight"/>.<see cref="IIndexable1D.Length"/></param>
        /// <param name="result">The vector whose entries will be overwritten with the result of the multiplication.</param>
        /// <param name="resultStart">The index of the first entry of <paramref name="result"/> that will be overwritten. 
        ///     Constraints: <paramref name="resultStart"/> + this.<see cref="NumRows"/> &lt;= 
        ///     <paramref name="result"/>.<see cref="IIndexable1D.Length"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the arguments do not satisfy the described 
        ///     constraints.</exception>
        public void MultiplyVectorSection(IVectorView vectorRight, int vectorStart, Vector result, int resultStart)
        {
            Preconditions.CheckMultiplicationDimensionsSection(this, vectorRight, vectorStart, result, resultStart);
            for (int i = 0; i < NumRows; ++i)
            {
                double dot = 0.0;
                int rowStart = rowOffsets[i]; //inclusive
                int rowEnd = rowOffsets[i + 1]; //exclusive
                for (int k = rowStart; k < rowEnd; ++k)
                {
                    dot += values[k] * vectorRight[vectorStart + colIndices[k]];
                }
                result[resultStart + i] = dot;
            }
        }

        /// <summary>
        /// See <see cref="IReducible.Reduce(double, ProcessEntry, ProcessZeros, Reduction.Finalize)"/>.
        /// </summary>
        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            int nnz = values.Length;
            for (int i = 0; i < nnz; ++i) aggregator = processEntry(values[i], aggregator);
            aggregator = processZeros(NumRows * NumColumns - nnz, aggregator);
            return finalize(aggregator);
        }

        /// <summary>
        /// See <see cref="IMatrixView.Scale(double)"/>.
        /// </summary>
        IMatrix IMatrixView.Scale(double scalar) => Scale(scalar);

        /// <summary>
        /// Performs the following operation for the non-zero entries (i, j), such that 0 &lt;= i &lt; <see cref="NumRows"/>, 
        /// 0 &lt;= j &lt; <see cref="NumColumns"/>: result[i, j] = <paramref name="scalar"/> * this[i, j].
        /// The resulting matrix is written to a new <see cref="CsrMatrix"/> and then returned.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this matrix.</param>
        public CsrMatrix Scale(double scalar)
        {
            int nnz = this.values.Length;
            double[] resultValues = new double[nnz];
            Array.Copy(this.values, resultValues, nnz); //TODO: perhaps I should also copy the indexers
            Blas.Dscal(nnz, scalar, resultValues, 0, 1);
            return new CsrMatrix(this.NumRows, this.NumColumns, resultValues, this.colIndices, this.rowOffsets);
        }

        /// <summary>
        /// See <see cref="IMatrix.ScaleIntoThis(double)"/>.
        /// </summary>
        public void ScaleIntoThis(double scalar) => Blas.Dscal(values.Length, scalar, values, 0, 1);

        /// <summary>
        /// See <see cref="IMatrix.SetEntryRespectingPattern(int, int, double)"/>.
        /// </summary>
        public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value)
        {
            int entryOffset = FindOffsetOf(rowIdx, colIdx);
            if (entryOffset == zeroEntryOffset) throw new SparsityPatternModifiedException(
                $"Cannot write to zero entry ({rowIdx}, {colIdx}).");
            else values[entryOffset] = value;
        }

        /// <summary>
        /// See <see cref="IMatrixView.Transpose"/>.
        /// </summary>
        public IMatrix Transpose() => TransposeToCSC(true);

        /// <summary>
        /// Creates a new <see cref="CscMatrix"/> instance, that is transpose to this: result[i, j] = this[j, i]. The 
        /// internal arrays can be copied or shared with this <see cref="CsrMatrix"/> instance.
        /// </summary>
        /// <param name="copyInternalArray">If true, the internal arrays that store the entries of this 
        ///     <see cref="CsrMatrix"/> instance will be copied and the new <see cref="CscMatrix"/> instance 
        ///     instance will have references to the copies, which is safer. If false, both the new matrix and this one will have  
        ///     references to the same internal arrays, which is faster.</param>
        public CscMatrix TransposeToCSC(bool copyInternalArrays)
        {
            if (copyInternalArrays)
            {
                double[] valuesCopy = new double[values.Length];
                Array.Copy(values, valuesCopy, values.Length);
                int[] colIndicesCopy = new int[colIndices.Length];
                Array.Copy(colIndices, colIndicesCopy, colIndices.Length);
                int[] rowOffsetsCopy = new int[rowOffsets.Length];
                Array.Copy(rowOffsets, rowOffsetsCopy, rowOffsets.Length);
                return CscMatrix.CreateFromArrays(NumColumns, NumRows, valuesCopy, colIndicesCopy, rowOffsetsCopy, false);
            }
            else return CscMatrix.CreateFromArrays(NumColumns, NumRows, values, colIndices, rowOffsets, false);
        }

        /// <summary>
        /// Creates a new <see cref="CsrMatrix"/> instance, that is transpose to this: result[i, j] = this[j, i].
        /// </summary>
        public CsrMatrix TransposeToCSR()
        {
            // Use C# port of the scipy method.
            // TODO: Perhaps it could be done faster by making extra assumptions. Otherwise use SparseBLAS
            int nnz = this.values.Length;
            var cscValues = new double[nnz];
            var cscRowIndices = new int[nnz];
            var cscColOffsets = new int[NumColumns + 1];

            SparseArrays.CsrToCsc(NumRows, NumColumns, this.rowOffsets, this.colIndices, this.values, 
                cscColOffsets, cscRowIndices, cscValues);

            return new CsrMatrix(NumColumns, NumRows, cscValues, cscRowIndices, cscColOffsets);
        }

        /// <summary>
        /// Return the index into values and colIndices arrays, if the (rowIdx, colIdx) entry is within the pattern. 
        /// Otherwise returns <see cref="zeroEntryOffset"/>.
        /// </summary>
        /// <param name="rowIdx"></param>
        /// <param name="colIdx"></param>
        private int FindOffsetOf(int rowIdx, int colIdx)
        {
            //TODO: if the column indices of each row are sorted, then this method is faster using binary search.
            Preconditions.CheckIndices(this, rowIdx, colIdx); //TODO: check indices?
            int rowStart = rowOffsets[rowIdx]; //inclusive
            int rowEnd = rowOffsets[rowIdx + 1]; //exclusive
            for (int k = rowStart; k < rowEnd; ++k)
            {
                if (colIndices[k] == colIdx) return k;
            }
            return zeroEntryOffset;
        }

        private bool HasSameIndexer(CsrMatrix other)
        {
            return (this.colIndices == other.colIndices) && (this.rowOffsets == other.rowOffsets);
        }
    }
}
